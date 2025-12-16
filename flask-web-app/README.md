# Flask Web App

This is an interactive web application built with Flask for visualizing and exploring biorepository and clinical datasets. The app provides dynamic filtering, bar charts, and donut plots using Altair, and allows users to drill down into patient-level details.

## Features
- **Interactive Bar Charts**: Visualize sample counts by cancer type or timepoint, with options to color by clinical or demographic variables.
- **Donut Plots**: Explore the distribution of material types (e.g., Plasma, Serum, PBMC) by patient count or tube count.
- **Barcodes and Patient Detail Pages**: View all samples and clinical metadata for individual patients.
- **Filtering**: Filter data by medical history, pre-treatment, irAE, cancer type, timepoint, and more.

## Folder Structure
- `app.py` — Main Flask application.
- `datasets/` — Contains CSV files with cleaned data for the app.
- `templates/` — HTML templates for rendering pages.
- `static/` — Static assets (CSS, JS, images).
- `requirements.txt` — Python dependencies.
- `Procfile` — For deployment (e.g., on Heroku).

## Setup 
1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
2. **Run the app**:
   ```bash
   python app.py
   ```
   The app will be available at `https://cponti.pythonanywhere.com/`.

## Requirements
- Python 3.7+
- Flask
- pandas
- altair
- matplotlib

Install all requirements with `pip install -r requirements.txt`.

## Deployment
- The app was deployed using pythoAnywhere: https://cponti.pythonanywhere.com/.


