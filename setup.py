from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

with open('requirements.txt') as f:
    lines = f.readlines()
install_requirements = [i.strip() for i in lines]


setup(name='isolines',
      version='0.1',
      description='street network isolines and isochrones',
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.6',
          'Topic :: Text Processing :: Linguistic',
      ],
      keywords='network isoline isochrone osm',
      url='http://github.com/mlichter2/isolines',
      author='mlichter2',
      author_email='mlichter@gmail.com',
      license='MIT',
      packages=['isolines'],
      install_requires=install_requirements,
      test_suite='nose.collector',
      tests_require=['nose', 'nose-cover3'],
      include_package_data=True,
      zip_safe=False)
