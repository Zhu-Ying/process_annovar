import setuptools

# with open("README.md", "r") as fh:
#     long_description = fh.read()

setuptools.setup(
    name="process_annovar",
    version="0.0.1",
    author="Ying Zhu",
    author_email="win19890412@163.com",
    # description="A tools of bam to image",
    # long_description=long_description,
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
