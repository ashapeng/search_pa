from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
import os

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/submit', methods=['POST'])
def submit():
    try:
        print("Received request")  # Debug print
        data = request.json
        if not data:
            return jsonify({"status": "error", "message": "No data received"}), 400
        
        name = data.get('name')
        school = data.get('school')
        program = data.get('program')
        citizenship = data.get('citizenship')
        
        # Here you can add logic to save or process the data
        print(f"Received: Name={name}, School={school}, Program={program}, Citizenship={citizenship}")
        
        return jsonify({"status": "success", "message": "Information received successfully!"})
    
    except Exception as e:
        print(f"Error processing request: {str(e)}")
        return jsonify({"status": "error", "message": str(e)}), 500

if __name__ == '__main__':
    # Use Waitress for production
    from waitress import serve
    print("Starting server on http://0.0.0.0:8000")
    serve(app, host="0.0.0.0", port=8000, url_scheme='http')