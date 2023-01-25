from setuptools import setup

setup(
    name='kintunGIP',
    version='0.9.0',
    packages=["kintunGIP"],
    scripts=["bin/kintunGIP"],
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
