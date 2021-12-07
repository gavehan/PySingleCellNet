from setuptools import setup

setup(name='pyParallelSingleCellNet',
      version='0.1',
      description='Determining cell identity from single cell RNA-Seq data',
      url='https://github.com/gavehan/pyParallelSingleCellNet/',
      author='Junhan Kim',
      author_email='junkim779@yonsei.ac.kr',
      license='MIT',
      packages=['pyParallelSingleCellNet'],
      install_requires=[
          'pandas',
          'numpy',
          'sklearn',
          'scanpy',
          'sklearn',
          'statsmodels',
          'scipy',
          'matplotlib',
          'seaborn',
          'umap-learn',
          'tqdm'
      ],
      zip_safe=False)
