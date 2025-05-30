from flask import Flask, render_template, request, jsonify
import json
import os
import requests
from bs4 import BeautifulSoup
import re
import time

app = Flask(__name__)

# Define the database file path
DB_FILE = 'mcgill_professors.json'

def load_database():
    print("\n=== Loading Database ===")
    if os.path.exists(DB_FILE):
        try:
            with open(DB_FILE, 'r') as f:
                data = json.load(f)
                if isinstance(data, list):
                    print(f"✓ Successfully loaded database with {len(data)} professors")
                    return data
                print("❌ Database format is incorrect")
                return []
        except json.JSONDecodeError:
            print("❌ Error decoding database file")
            return []
    print("❌ Database file does not exist")
    return []

def save_database(db):
    print("\n=== Saving Database ===")
    try:
        with open(DB_FILE, 'w') as f:
            json.dump(db, f, indent=4)
        print(f"✓ Successfully saved database with {len(db)} professors")
    except Exception as e:
        print(f"❌ Error saving database: {str(e)}")

def scrape_mcgill_biology():
    try:
        print("\n=== Starting McGill Biology Scraping Process ===")
        
        # Step 1: Get the search term
        search_term = request.get_json().get('search_term', '').strip()
        if not search_term:
            print("❌ No search term provided")
            return []
        
        print(f"Searching for professor: {search_term}")
        
        # Setup session with browser-like headers
        session = requests.Session()
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Accept-Language': 'en-US,en;q=0.5',
            'Referer': 'https://www.mcgill.ca/biology/people',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1',
            'Cache-Control': 'max-age=0',
            'Host': 'www.mcgill.ca'
        }
        session.headers.update(headers)

        # Function to clean professor name
        def clean_name(name):
            name = name.strip()
            name = re.sub(r'\s*\(.*?\)', '', name)  # Remove parentheses and content
            name = re.sub(r'\s*Ph\.?D\.?', '', name, flags=re.IGNORECASE)  # Remove PhD
            name = re.sub(r'^(Dr\.|Professor|Prof\.)\s+', '', name, flags=re.IGNORECASE)  # Remove titles
            return name

        # Function to check if names match
        def names_match(search_name, candidate_name):
            search_parts = search_name.lower().split()
            candidate_parts = candidate_name.lower().split()
            return all(part in ' '.join(candidate_parts) for part in search_parts)

        # Create URL-friendly version of the name
        url_name = search_term.lower().replace(' ', '-')
        
        # List of possible URL patterns
        urls_to_try = [
            f"https://www.mcgill.ca/biology/{url_name}",
            f"https://www.mcgill.ca/biology/people/{url_name}",
            f"https://www.mcgill.ca/biology/staff/{url_name}",
            f"https://www.mcgill.ca/biology/faculty/{url_name}"
        ]

        professors = []
        
        # Try each URL pattern
        for url in urls_to_try:
            try:
                print(f"\nTrying URL: {url}")
                response = session.get(url, allow_redirects=True)
                
                if response.status_code == 200:
                    print("✓ URL access successful")
                    
                    # Save response for debugging
                    with open(f'response_{len(professors)}.html', 'w', encoding='utf-8') as f:
                        f.write(response.text)
                    
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Try different methods to find the professor name
                    name_found = False
                    professor_name = None
                    
                    # Method 1: Page title
                    title = soup.find('title')
                    if title:
                        title_text = title.text.strip()
                        if 'Biology' in title_text and not any(x in title_text.lower() for x in ['not found', 'error', 'rejected']):
                            professor_name = title_text.split('|')[0].strip()
                            name_found = True
                    
                    # Method 2: Main heading
                    if not name_found:
                        heading = soup.find(['h1', 'h2'], class_=['page-header', 'title'])
                        if heading:
                            professor_name = heading.text.strip()
                            name_found = True
                    
                    if professor_name:
                        professor_name = clean_name(professor_name)
                        if names_match(search_term, professor_name):
                            print(f"✓ Found matching professor: {professor_name}")
                            professors.append({
                                'name': professor_name,
                                'url': url,
                                'department': 'Biology',
                                'university': 'McGill University'
                            })
                            print(f"✓ Added professor with URL: {url}")
                            break  # Found the professor, no need to try other URLs
                
            except Exception as e:
                print(f"❌ Error accessing URL: {str(e)}")
                continue

        if not professors:
            # If no professors found, try a simpler approach
            print("\nTrying direct URL...")
            direct_url = f"https://www.mcgill.ca/biology/{url_name}"
            professors.append({
                'name': search_term,
                'url': direct_url,
                'department': 'Biology',
                'university': 'McGill University'
            })
            print(f"Added professor using direct URL: {direct_url}")

        print(f"\n=== Search Complete ===")
        print(f"🎯 Total professors found: {len(professors)}")
        return professors
            
    except Exception as e:
        print(f"❌ Error in main scraping process: {str(e)}")
        import traceback
        print(f"Stack trace: {traceback.format_exc()}")
        return []

# Load the database when the application starts
professors_db = load_database()

@app.route('/')
def index():
    return render_template('network_prof.html')

@app.route('/search_professors', methods=['POST'])
def search_professors():
    print("\n=== Processing Search Request ===")
    data = request.get_json()
    search_term = data.get('search_term', '').strip()
    department = data.get('department', 'Biology').lower()
    
    if not search_term:
        return jsonify({
            'professors': [],
            'total': 0,
            'error': 'Please enter a search term'
        })
    
    print(f"Search term: {search_term}")
    print(f"Department: {department}")
    
    # Always scrape for professors
    print("Starting scraping process...")
    scraped_professors = scrape_mcgill_biology()
    
    # Search through scraped professors with more flexible matching
    results = []
    search_parts = search_term.lower().split()
    
    for prof in scraped_professors:
        prof_name = prof['name'].lower()
        # Check if all parts of the search term are in the professor's name
        if all(part in prof_name for part in search_parts):
            results.append(prof)
            print(f"✓ Found matching professor: {prof['name']}")
    
    print(f"Total results found: {len(results)}")
    
    if not results:
        # Save the response for debugging
        print("No results found. Saving debug information...")
        with open('debug_search.txt', 'a') as f:
            f.write(f"\n\nSearch attempt at {time.strftime('%Y-%m-%d %H:%M:%S')}:\n")
            f.write(f"Search term: {search_term}\n")
            f.write(f"Department: {department}\n")
            f.write(f"Total professors scraped: {len(scraped_professors)}\n")
            f.write("Scraped professors:\n")
            for prof in scraped_professors:
                f.write(f"- {prof['name']} ({prof['url']})\n")
    
    return jsonify({
        'professors': results,
        'total': len(results)
    })

@app.route('/get_professor_details', methods=['POST'])
def get_professor_details():
    data = request.get_json()
    name = data.get('name')
    
    if not name:
        return jsonify({'error': 'Name is required'}), 400
    
    # Find the professor
    professor = None
    for prof in professors_db:
        if prof['name'].lower() == name.lower():
            professor = prof
            break
    
    if not professor:
        return jsonify({'error': 'Professor not found'}), 404
    
    # Find related professors (same department)
    connections = []
    for prof in professors_db:
        if (prof['name'].lower() != name.lower() and 
            prof['department'].lower() == professor['department'].lower()):
            connections.append(prof)
    
    return jsonify({
        'professor': professor,
        'connections': connections
    })

if __name__ == '__main__':
    app.run(debug=True) 