import setuptools

with open("README.org", r) as fh:
    long_description = fh.read()

setuptools.setup(
    name="amandin-Zain-Jabbar",
    version="0.0.1",
    author="Zain Jabbar",
    author_email="zaijab2000@gmail.com",
    description="Data Analysis Package for Manoa",
    long_description=long_description,
    long_description_content_type="text/plain-text",
    url="https://github.com/zaijab/amandin",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
