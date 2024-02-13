from setuptools import setup, find_packages 



with open("README.md", "r") as f:
    long_description = f.read()
    
    
    
setup(
    name="interface_builder",
    version='0.0.1',
    author='',
    author_email='',
    description="build silica-water interface",
    license='BSD',
    long_description=long_description, 
    keywords='silica,molecular dynamics',
    url="", 
    packages=find_packages(),
    install_requires=[
        "numpy",
        "molecular-builder"
    ],
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "Development Status :: 2 - Pre-Alpha",
        "Operating System :: Unix",
    ]
)
