import requests
import webbrowser
import argparse
import os
from bs4 import BeautifulSoup
import re
import PyPDF2

def fetch_mcgill_handbook(filename:str):
    # URL of the McGill Biology Department's graduate page
    url = 'https://www.mcgill.ca/biology/graduate-0/current-graduate-students-0'

    try:
        # Send a GET request to the page
        response = requests.get(url)
        response.raise_for_status()  # Check for request errors

        # Parse the page content
        soup = BeautifulSoup(response.text, 'html.parser')

        # Define a vague search pattern for "handbook"
        search_pattern = re.compile(r'handbook', re.IGNORECASE)

        # Find all paragraphs and list items containing the word "handbook"
        relevant_sections = soup.find_all(lambda tag: tag.name in ["p", "li"] and search_pattern.search(tag.get_text()))

        # Extract links from those sections
        matching_links = []
        for section in relevant_sections:
            links = section.find_all('a', href=True)
            for link in links:
                href = link['href']
                # Convert relative URLs to absolute
                if not href.startswith('http'):
                    href = 'https://www.mcgill.ca' + href
                matching_links.append((section.get_text(strip=True), href))

        # Display results， if no results, print "No handbook-related links found on the page." and return None and exit the function
        if matching_links:
            print("Handbook-related links found")
            print(f"Found {len(matching_links)} matching links")
        
            # check if the link is a PDF, and store PDFs in a list
            pdf_links = []
            for text, href in matching_links:
                if href.endswith('.pdf'):
                    pdf_url = href
                    pdf_links.append(pdf_url)
            
            # if multiple PDFs are found, prompt the user to choose one or multiple PDFs to download
            if len(pdf_links) > 1:
                print("Multiple PDFs found:")
                for i, link in enumerate(pdf_links):
                    print(f"{i+1}. {link}")
                choices = int(input("Enter all the numbers of the PDF you want to download: "))
                pdf_urls = [pdf_links[i-1] for i in choices]
            else:
                pdf_urls = [pdf_links[0]]

            # check if the filename has the .pdf extension, if yes, remove it
            if filename.endswith(".pdf"):
                filename = filename[:-4]

            # Download the PDF and save
            downloaded_files = []
            for i, pdf_url in enumerate(pdf_urls):
                pdf_filename = f"{filename}_{i+1}.pdf"
                
                with open(pdf_filename, 'wb') as f:
                    f.write(requests.get(pdf_url).content)
                print(f"PDF saved to {pdf_filename}")
                downloaded_files.append(pdf_filename)

                # Open the PDF in the default PDF viewer
                webbrowser.open(pdf_filename)
            
            return downloaded_files

        else:
            print("No handbook-related links found on the page.")
            return None
    
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        return None

# define a function to search for a specific topic in the handbook
def search_handbook(topic:str, pdf_files:list):
    # read the handbook
    for pdf_file in pdf_files:
        try:
            with open(pdf_file, "rb") as f:
                pdf_reader = PyPDF2.PdfReader(f)
                # search for the topic in the handbook
                for page in pdf_reader.pages:
                    text = page.extract_text()
                    if topic.lower() in text.lower():
                        print(f"\nFound in {pdf_file}:")
                        print(text)
        except Exception as e:
            print(f"Error processing {pdf_file}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Fetch the McGill Biology Department's graduate handbook")
    parser.add_argument("-o", "--output", help="Output directory and main file name (default: the current py directory)")
    
    args = parser.parse_args()

    # Set default output directory and main file name if not specified
    if not args.output:
        # prompt the user to enter a directory to save the PDF files
        args.output = input("Enter a directory to save the PDF files: ")
    else:
        args.output = os.path.join(args.output)

    pdf_files = fetch_mcgill_handbook(args.output)

    # Only proceed with search if we have PDF files
    if pdf_files:
        # search for a specific topic in the handbook
        # prompt the user to enter a topic to search for
        topic = input("Enter a topic to search for: ")
        search_handbook(topic, pdf_files)
    else:
        print("No PDF files were downloaded. Cannot perform search.")



if __name__ == "__main__":
    main()
