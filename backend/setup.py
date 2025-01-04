from setuptools import setup, find_packages

setup(
    name="chemical_calculator",  # Replace with your project name
    version="0.1.0",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "chemical_calculator=app:app",  # This will point to the `app` object inside `app.py`
        ],
    },
    install_requires=[
    "Flask",  # If you're building a web application
    "numpy",  # For numerical computations (if used in your project)
    "pandas",  # Data handling (if applicable)
    "scipy",  # For scientific computations (if needed)
    "matplotlib",  # For plotting or visualizations (if applicable)
    "requests",  # If you're making HTTP requests
    "flask_sqlalchemy",  # If you're using a database in Flask
    "pytest",  # For testing (if you're writing tests)
    "flask-cors",  # If you need Cross-Origin Resource Sharing support in Flask
    "openpyxl",  # If you're working with Excel files (if applicable)
    "pillow",  # For image processing (if used in your project)
    "jsonschema",  # If you're working with JSON data and schemas
    "pyyaml",  # For YAML handling (if needed)
    "flask-wtf",  # If you're using Flask forms
    "Flask-Login",  # If your project has user authentication
    "flask-migrate",  # If you're doing database migrations in Flask
    # Add more dependencies as necessary
],

    description="A backend for the chemical calculator project",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/your-repo-link",  # Optional
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
)
