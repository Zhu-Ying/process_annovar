import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="annovar_tools",
    version="0.0.7",
    author="Ying Zhu",
    author_email="win19890412@163.com",
    description="tools for ANNOVAR",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'pysam',
    ],
    scripts=['annovar_tools.py']
)
