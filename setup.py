from setuptools import find_packages, setup


setup(
    name="dorfannotator",
    version="1.0.0",
    description="Functional annotation of variants affecting downstream open reading frames in 3′UTRs.",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    license="GPL-3.0-only",
    author="Sihan Liu",
    python_requires=">=3.10",
    packages=find_packages(exclude=("tests", "tests.*")),
    install_requires=[],
    entry_points={
        "console_scripts": [
            "dORFannotator=dorfannotator.cli:main",
        ],
    },
)
