from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()
    

setup(name='pfpostproc',
      version='0.0.1',
      description='parflow-clm postprocessing tools',
      long_description=readme(),
      long_description_content_type='text/markdown',
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3',
          'Operating System :: OS Independent'
      ],
      url='https://github.com/luketelfer/pfpostproc',
      author='Luke Telfer',
      author_email='luketelfer@boisestate.edu',
      keywords='parflow clm postproc',
      license='MIT',
      packages=['pfpostproc'],
      install_requires=[],
      include_package_data=True,
      zip_safe=False)