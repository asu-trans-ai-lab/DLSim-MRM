import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

try:
    # if have requirements.txt file inside the folder
    with open("requirements.txt", "r", encoding="utf-8") as f:
        modules_needed = [i.strip() for i in fh.readlines()]
except Exception:
    modules_needed = []

setuptools.setup(
    name="DLSim",
    version="0.2.11",
    author="Dr.Xuesong (Simon) Zhou, Dr.Cafer Avci, Xiangyong Luo",
    author_email="xzhou74@asu.edu",
    License="Apache Software License (http://www.apache.org/licenses/LICENSE-2.0)",
    description="DLSim is an open-source, cross-platform, lightweight, and fast Python traffic assignment tool adopted and modified from ASU TransAI Lab",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/asu-trans-ai-lab/DLSim",
    install_requires=modules_needed,
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: Implementation :: PyPy",
    ],

    packages=setuptools.find_packages(where='DLSim'),
    package_dir={'': 'DLSim'},
    include_package_data=True,

    package_data={'': ['*.txt', '*.xls', '*.xlsx', '*.csv', '*.png', "*.dll", "*.so", "*.dylib"],
                  "pydtalite_bin": ["*.dll", "*.so", "*.dylib"]},
)