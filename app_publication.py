from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
import requests
from datetime import datetime
import time
from scholarly import scholarly
from Bio import Entrez
import json
import re

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes

# Configure Entrez email (required for PubMed API)
Entrez.email = "your_email@example.com"  # Replace with your email

RESEARCH_FIELDS = {
    'computer_science': {
        'keywords': ['computer science', 'artificial intelligence', 'machine learning', 'data mining', 'software', 
                    'algorithm', 'computing', 'programming', 'database', 'network', 'cybersecurity', 
                    'information system', 'cloud computing', 'distributed system'],
        'journals': ['IEEE', 'ACM', 'Computer', 'Computing', 'Informatics', 'Software']
    },
    'biology': {
        'keywords': ['biology', 'molecular', 'cell', 'genetics', 'genome', 'protein', 'rna', 'dna',
                    'organism', 'species', 'evolution', 'enzyme', 'biological'],
        'journals': ['Cell', 'Nature', 'Science', 'Molecular Biology', 'Genetics', 'Development']
    },
    'chemistry': {
        'keywords': ['chemistry', 'chemical', 'molecule', 'synthesis', 'reaction', 'compound', 'polymer',
                    'catalyst', 'organic', 'inorganic', 'biochemistry'],
        'journals': ['Chemical', 'Chemistry', 'Journal of the American Chemical Society', 'Angewandte Chemie']
    },
    'physics': {
        'keywords': ['physics', 'quantum', 'particle', 'mechanics', 'relativity', 'electromagnetic',
                    'nuclear', 'optics', 'thermodynamics', 'condensed matter'],
        'journals': ['Physical Review', 'Physics', 'Nature Physics', 'Journal of Physics']
    },
    'mathematics': {
        'keywords': ['mathematics', 'mathematical', 'algebra', 'geometry', 'topology', 'analysis',
                    'probability', 'statistics', 'theorem', 'numerical'],
        'journals': ['Mathematical', 'Mathematics', 'Journal of Algebra', 'Topology', 'Statistics']
    },
    'engineering': {
        'keywords': ['engineering', 'mechanical', 'electrical', 'civil', 'chemical engineering',
                    'aerospace', 'robotics', 'control system', 'materials science'],
        'journals': ['Engineering', 'Journal of Engineering', 'IEEE Transactions']
    }
}

def normalize_institution(institution):
    """Normalize institution name for better matching"""
    replacements = {
        "university of": "univ",
        "university": "univ",
        "of": "",
    }
    
    name = institution.lower()
    for old, new in replacements.items():
        name = name.replace(old, new)
    return ' '.join(word for word in name.split() if word)

def search_google_scholar(professor_name, institution, year_start, year_end):
    """Search for publications using Google Scholar"""
    try:
        print(f"Searching Google Scholar for: {professor_name}")
        
        # Search for the author
        search_query = scholarly.search_author(f"{professor_name} {institution}")
        author = next(search_query)
        
        # Fill in author data
        author = scholarly.fill(author)
        
        # Get publications
        publications = []
        for pub in author['publications']:
            try:
                pub_complete = scholarly.fill(pub)
                year = pub_complete.get('bib', {}).get('pub_year')
                
                if year and year_start <= int(year) <= year_end:
                    publications.append({
                        'title': pub_complete['bib'].get('title', 'No title'),
                        'year': int(year),
                        'url': f"https://scholar.google.com/citations?view_op=view_citation&citation_for_view={pub_complete['author_pub_id']}",
                        'venue': pub_complete['bib'].get('venue', 'No venue available'),
                        'abstract': pub_complete.get('bib', {}).get('abstract', 'No abstract available'),
                        'source': 'Google Scholar'
                    })
            except Exception as e:
                print(f"Error processing publication: {str(e)}")
                continue
        
        return publications
    except Exception as e:
        print(f"Google Scholar Error: {str(e)}")
        return []

def search_pubmed(professor_name, institution, year_start, year_end, email=None, author_position='any'):
    """Search for publications using PubMed"""
    try:
        if not email:
            print("No email provided for PubMed search")
            return []
            
        print(f"Searching PubMed for: {professor_name}")
        Entrez.email = email  # Set email for this search
        
        # Search for author's papers
        search_term = f"{professor_name}[Author] AND (\"{institution}\"[Affiliation])"
        if author_position == 'first':
            search_term = f"{professor_name}[1st Author] AND (\"{institution}\"[Affiliation])"
        
        handle = Entrez.esearch(db="pubmed", term=search_term, 
                              mindate=str(year_start), maxdate=str(year_end),
                              retmax=100)
        record = Entrez.read(handle)
        handle.close()
        
        if not record['IdList']:
            return []
            
        # Fetch paper details
        handle = Entrez.efetch(db="pubmed", id=record['IdList'], rettype="medline", retmode="text")
        records = Entrez.parse(handle)
        
        publications = []
        for record in records:
            try:
                year = int(record.get('DP', '0')[:4])
                if not (year_start <= year <= year_end):
                    continue

                # Get authors list
                authors = record.get('AU', [])
                if not authors:
                    continue

                # Check author positions
                is_first_author = authors[0].lower() == professor_name.lower()
                # For PubMed, check if professor is in the last position (typical for corresponding author)
                is_corresponding_author = authors[-1].lower() == professor_name.lower()

                # Apply author position filter
                if author_position == 'first' and not is_first_author:
                    continue
                if author_position == 'corresponding' and not is_corresponding_author:
                    continue
                if author_position == 'both' and not (is_first_author and is_corresponding_author):
                    continue

                publications.append({
                    'title': record.get('TI', 'No title'),
                    'year': year,
                    'url': f"https://pubmed.ncbi.nlm.nih.gov/{record.get('PMID', '')}",
                    'venue': record.get('TA', 'No venue available'),
                    'abstract': record.get('AB', 'No abstract available'),
                    'source': 'PubMed',
                    'is_first_author': is_first_author,
                    'is_corresponding_author': is_corresponding_author
                })
            except Exception as e:
                print(f"Error processing PubMed record: {str(e)}")
                continue
                
        return publications
    except Exception as e:
        print(f"PubMed Error: {str(e)}")
        return []

def search_semantic_scholar(professor_name, institution, year_start, year_end, author_position='any'):
    """Search for publications using Semantic Scholar"""
    base_url = "https://api.semanticscholar.org/graph/v1"
    
    try:
        # Search for the author
        search_url = f"{base_url}/author/search"
        params = {
            "query": professor_name,
            "fields": "authorId,name,affiliations,url,paperCount"
        }
        
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        author_data = response.json()
        
        if not author_data.get("data"):
            return []
            
        # Find best match
        best_match = None
        for author in author_data["data"]:
            if author.get('name', '').lower() == professor_name.lower():
                best_match = author
                break
        
        if not best_match:
            return []
            
        # Get author's papers
        papers_url = f"{base_url}/author/{best_match['authorId']}/papers"
        params = {
            "fields": "title,year,url,abstract,venue,authors",
            "limit": 100
        }
        
        papers_response = requests.get(papers_url, params=params)
        papers_response.raise_for_status()
        papers_data = papers_response.json()
        
        publications = []
        for paper in papers_data.get("data", []):
            year = paper.get("year")
            if not year or not (year_start <= year <= year_end):
                continue

            authors = paper.get("authors", [])
            if not authors:
                continue

            # Check if the professor is first author or corresponding author
            is_first_author = authors[0].get('name', '').lower() == professor_name.lower()
            # For Semantic Scholar, we'll assume the last author is the corresponding author
            is_corresponding_author = authors[-1].get('name', '').lower() == professor_name.lower()

            # Apply author position filter
            if author_position == 'first' and not is_first_author:
                continue
            if author_position == 'corresponding' and not is_corresponding_author:
                continue
            if author_position == 'both' and not (is_first_author and is_corresponding_author):
                continue

            publications.append({
                "title": paper.get("title", "No title"),
                "year": year,
                "url": paper.get("url", "No URL available"),
                "venue": paper.get("venue", "No venue available"),
                "abstract": paper.get("abstract", "No abstract available"),
                "source": "Semantic Scholar",
                "is_first_author": is_first_author,
                "is_corresponding_author": is_corresponding_author
            })
        
        return publications
    except Exception as e:
        print(f"Semantic Scholar Error: {str(e)}")
        return []

def merge_and_deduplicate_papers(papers_list):
    """Merge papers from different sources and remove duplicates"""
    seen_titles = set()
    merged_papers = []
    
    for paper in papers_list:
        # Normalize title for comparison
        normalized_title = re.sub(r'[^a-zA-Z0-9\s]', '', paper['title'].lower())
        
        if normalized_title not in seen_titles:
            seen_titles.add(normalized_title)
            merged_papers.append(paper)
    
    return merged_papers

def filter_papers_by_field(papers, target_field):
    """Filter papers based on research field"""
    if not target_field or target_field not in RESEARCH_FIELDS:
        return papers
    
    field_info = RESEARCH_FIELDS[target_field]
    keywords = [kw.lower() for kw in field_info['keywords']]
    journals = [j.lower() for j in field_info['journals']]
    
    filtered_papers = []
    for paper in papers:
        # Check title, abstract, and venue for field relevance
        title = paper.get('title', '').lower()
        abstract = paper.get('abstract', '').lower()
        venue = paper.get('venue', '').lower()
        
        # Check if any keyword is present in title or abstract
        keyword_match = any(kw in title or kw in abstract for kw in keywords)
        # Check if any journal keyword is present in venue
        journal_match = any(j in venue for j in journals)
        
        if keyword_match or journal_match:
            filtered_papers.append(paper)
    
    return filtered_papers

def search_publications(professor_name, institution, year_start=None, year_end=None, email=None, 
                       first_author_only=False, corresponding_author_only=False, both_positions=False,
                       research_field=None):
    """Search for publications across all available sources"""
    # Set default year range to last 5 years
    if year_start is None or year_end is None:
        current_year = datetime.now().year
        year_end = current_year
        year_start = current_year - 5
    
    # Determine author position filter
    author_position = 'any'
    if both_positions:
        author_position = 'both'
    elif first_author_only:
        author_position = 'first'
    elif corresponding_author_only:
        author_position = 'corresponding'
    
    all_papers = []
    
    # Start with PubMed if email is provided
    if email:
        print(f"Searching PubMed for: {professor_name}")
        pubmed_papers = search_pubmed(professor_name, institution, year_start, year_end, email, author_position)
        if research_field:
            pubmed_papers = filter_papers_by_field(pubmed_papers, research_field)
        all_papers.extend(pubmed_papers)
        print(f"Found {len(pubmed_papers)} papers from PubMed")
        
        # If we get enough results from PubMed (e.g., more than 10 papers), just use those
        if len(pubmed_papers) > 10:
            print("Sufficient results from PubMed - skipping other sources")
            return {"papers": pubmed_papers}
    else:
        print("Skipping PubMed search - no email provided")
    
    # Try Semantic Scholar next
    print(f"Searching Semantic Scholar for: {professor_name}")
    semantic_papers = search_semantic_scholar(professor_name, institution, year_start, year_end, author_position)
    if research_field:
        semantic_papers = filter_papers_by_field(semantic_papers, research_field)
    all_papers.extend(semantic_papers)
    print(f"Added {len(semantic_papers)} papers from Semantic Scholar")
    
    # Finally try Google Scholar if we still need more results
    if len(all_papers) < 10:
        print(f"Searching Google Scholar for: {professor_name}")
        scholar_papers = search_google_scholar(professor_name, institution, year_start, year_end)
        if research_field and author_position == 'any':
            scholar_papers = filter_papers_by_field(scholar_papers, research_field)
            all_papers.extend(scholar_papers)
            print(f"Added {len(scholar_papers)} papers from Google Scholar")
        else:
            print("Skipping Google Scholar results due to author position filtering")
    
    # Merge and deduplicate
    merged_papers = merge_and_deduplicate_papers(all_papers)
    
    # Sort by year (newest first)
    merged_papers.sort(key=lambda x: x['year'], reverse=True)
    
    if not merged_papers:
        return {"error": f"No papers found for {professor_name} in {research_field if research_field else 'any field'} between {year_start} and {year_end}"}
    
    print(f"Final total: {len(merged_papers)} unique papers")
    return {"papers": merged_papers}

@app.route('/')
def index():
    return render_template('publication.html')

@app.route('/fields', methods=['GET'])
def get_fields():
    """Return available research fields"""
    return jsonify({
        'fields': [
            {'value': field, 'label': field.replace('_', ' ').title()} 
            for field in RESEARCH_FIELDS.keys()
        ]
    })

@app.route('/search', methods=['POST'])
def search():
    try:
        data = request.json
        professor_name = data.get('professor_name')
        institution = data.get('institution')
        email = data.get('email')
        year_start = data.get('year_start')
        year_end = data.get('year_end')
        first_author_only = data.get('first_author_only', False)
        corresponding_author_only = data.get('corresponding_author_only', False)
        both_positions = data.get('both_positions', False)
        research_field = data.get('research_field')
        
        # Validate inputs
        if not professor_name or not institution:
            return jsonify({"error": "Please provide both professor name and institution"})

        # Convert year strings to integers if provided
        if year_start:
            try:
                year_start = int(year_start)
            except ValueError:
                return jsonify({"error": "Invalid start year format"})

        if year_end:
            try:
                year_end = int(year_end)
            except ValueError:
                return jsonify({"error": "Invalid end year format"})

        # Validate years if provided
        current_year = datetime.now().year
        if year_start and year_end:
            if year_start > year_end:
                return jsonify({"error": "Start year cannot be greater than end year"})
            if year_start < 1900 or year_end > current_year:
                return jsonify({"error": f"Please enter years between 1900 and {current_year}"})
        elif year_start is None and year_end is None:
            # Use default 5-year range if no years provided
            year_end = current_year
            year_start = current_year - 5
        else:
            return jsonify({"error": "Please provide both start and end years or leave both empty"})
        
        results = search_publications(
            professor_name, 
            institution, 
            year_start, 
            year_end, 
            email,
            first_author_only,
            corresponding_author_only,
            both_positions,
            research_field
        )
        return jsonify(results)
        
    except Exception as e:
        return jsonify({"error": f"An unexpected error occurred: {str(e)}"})

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5000, debug=True)