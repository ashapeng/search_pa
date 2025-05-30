from flask import Flask, render_template, request, jsonify
import json
import os

app = Flask(__name__)

# Define the database file path
DB_FILE = 'people_database.json'

def load_database():
    if os.path.exists(DB_FILE):
        try:
            with open(DB_FILE, 'r') as f:
                data = json.load(f)
                # Ensure we have a list of dictionaries
                if isinstance(data, list):
                    return data
                return []
        except json.JSONDecodeError:
            return []
    return []

def save_database(db):
    with open(DB_FILE, 'w') as f:
        json.dump(db, f, indent=4)

# Load the database when the application starts
people_db = load_database()

@app.route('/')
def index():
    return render_template('network.html')

@app.route('/search_network', methods=['POST'])
def search_network():
    data = request.get_json()
    name = data.get('name')
    company = data.get('company')
    role = data.get('role')
    supervisor_name = data.get('supervisor_name')
    lab_website = data.get('lab_website')

    if not name or not company or not role:
        return jsonify({'error': 'Name, company, and role are required'}), 400

    # Search for the person in the database
    person = None
    for p in people_db:
        if isinstance(p, dict) and p.get('name', '').lower() == name.lower() and p.get('company', '').lower() == company.lower():
            person = p
            break

    # If person not found, create a new entry
    is_new = False
    if not person:
        is_new = True
        person = {
            'name': name,
            'company': company,
            'role': role,
            'supervisor_name': supervisor_name,
            'lab_website': lab_website,
            'connections': []
        }
        people_db.append(person)
        save_database(people_db)

    # Find connections (people in the same company)
    connections = []
    for p in people_db:
        if isinstance(p, dict) and p.get('name', '').lower() != name.lower() and p.get('company', '').lower() == company.lower():
            connections.append(p)

    return jsonify({
        'person': person,
        'connections': connections,
        'is_new': is_new
    })

@app.route('/update_person', methods=['POST'])
def update_person():
    data = request.get_json()
    name = data.get('name')
    supervisor_name = data.get('supervisor_name')
    lab_website = data.get('lab_website')

    if not name:
        return jsonify({'error': 'Name is required'}), 400

    # Find and update the person
    person = None
    for p in people_db:
        if isinstance(p, dict) and p.get('name', '').lower() == name.lower():
            person = p
            if supervisor_name:
                p['supervisor_name'] = supervisor_name
            if lab_website:
                p['lab_website'] = lab_website
            break

    if not person:
        return jsonify({'error': 'Person not found'}), 404

    save_database(people_db)
    return jsonify({'message': 'Person updated successfully', 'person': person})

if __name__ == '__main__':
    app.run(debug=True)