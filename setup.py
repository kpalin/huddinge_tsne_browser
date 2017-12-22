from setuptools import setup

requirements = [
    # package requirements go here
]

setup(
    name='huddinge_tsne_browser',
    version='0.1.0',
    description="Tool to browse sequence kmers laid out in 2D with TSNE approximating Huddinge distance",
    author="Kimmo Palin",
    author_email='kimmo.palin@helsinki.fi',
    url='https://github.com/kpalin/huddinge_tsne_browser',
    packages=['huddinge_tsne_browser'],
    entry_points={
        'console_scripts': [
            'huddinge_tsne_browser=huddinge_tsne_browser.cli:cli'
        ]
    },
    install_requires=requirements,
    keywords='huddinge_tsne_browser',
    classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6',
    ]
)
