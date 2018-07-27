from setuptools import setup

requirements = [
    "setuptools"
    # package requirements go here
]

setup(
    name='huddinge_tsne_browser',
    version='0.2.1',
    description="Tool to browse sequence kmers laid out in 2D with TSNE approximating Huddinge distance",
    author="Kimmo Palin",
    author_email='kimmo.palin@helsinki.fi',
    url='https://github.com/kpalin/huddinge_tsne_browser',
    packages=['huddinge_tsne_browser'],
    entry_points={
        'console_scripts':
        ['huddinge_tsne_browser=huddinge_tsne_browser.main:main']
    },
    #package_data={"huddinge_tsne_browser": ["resources/*.tsne"]},
    install_requires=requirements,
    keywords='huddinge_tsne_browser',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ])
