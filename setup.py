from setuptools import setup, find_packages

setup(
    name='kintunGIP',
    version='0.9.0',
    packages=find_packages(include=["kintunGIP", "kintunGIP.*"]),
    scripts=["kintunGIP/kintunGIP.py"],
    install_requires=[
        'biopython>=1.80',
        'glob2>=0.7',
    ],
    url='https://github.com/GMI-Lab/Kintun-GIP',
    license='GNU GENERAL PUBLIC LICENSE',
    author='Camilo Berrios-Pasten',
    author_email='camilo.berrios.p@ug.uchile.cl',
    description='Bacterial Genomic Island Predictor'
)
