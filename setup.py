from setuptools import setup, find_packages

setup(
    name='solar_alignment_py',
    version='1.0',
    description='Library to Align Solar Observations ALMA, IRIS, SDO',
    author='F.J Ordonez A. and J.C Guevara G',
    url='https://github.com/JavierOrdonezA/Alignment-of-IRIS-ALMA-and-SDO-images.git', 
    packages=find_packages(),
    author_email='fordonezaraujo@gmail.com -- juancamilo.guevaragomez@gmail.com',
    long_description='https://arxiv.org/abs/2404.04401 ', #Article 
    description_thesis= 'https://repositorio.unal.edu.co/handle/unal/85838', # Second chapter of my master thesis 
    python_requires='>=3.6',
)
