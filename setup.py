from setuptools import setup, find_packages

setup(
    name='scConGraph',
    version='0.1.12',
    description="a package for sc-RNA",
    long_description=open('README.md').read(),
    include_package_data=True,
    author='XinQi Li',
    author_email='lxq19@mails.tsinghua.edu.cn',
    maintainer='XinQi Li',
    maintainer_email='lxq19@mails.tsinghua.edu.cn',
    license='MIT License',
    url='',
    packages=find_packages(include=['scConGraph', 'scConGraph.*']),
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
    ],
    python_requires='>=3.7',
    install_requires=['network','matplotlib','numpy','pandas','seaborn','scipy','scanpy'],
)