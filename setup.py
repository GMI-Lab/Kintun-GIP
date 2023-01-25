from setuptools import setup, find_packages

setup(
    name='kintunGIP',
    version='0.9.0',
    packages=find_packages(include=["kintunGIP", "kintunGIP.*"]),
    install_requires=[
        'aragorn>=1.2.41',
        'barrnap>=0.9',
        'biopython>=1.80',
        'glob2>=0.7',
        'mmseqs2>=14.7e284',
        'python>=3.10',
        'sibeliaz>=1.2.5'
    ],
    url='https://github.com/GMI-Lab/Kintun-GIP',
    license='GNU GENERAL PUBLIC LICENSE',
    author='Camilo Berrios-Pasten',
    author_email='camilo.berrios.p@ug.uchile.cl',
    description='Bacterial Genomic Island Predictor'
)
