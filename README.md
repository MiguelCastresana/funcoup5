# FunCoup5 Network Analysis

This repository contains scripts and tools for conducting network analysis using **FunCoup5**, a comprehensive method for inferring functional coupling between genes and proteins.

## Overview

**FunCoup5** is a powerful tool for identifying functional associations between genes and proteins based on various data sources. It integrates diverse types of evidence, such as co-expression, protein-protein interactions, and phylogenetic profiles, to predict functional relationships within biological networks. The method has been widely used in bioinformatics and computational biology to explore the underlying functional organization of genomes.

Access the database here: [FunCoup5 website](https://funcoup.org/search/)
For a detailed explanation of **FunCoup5** and its applications, please refer to the [FunCoup5 paper](https://pubmed.ncbi.nlm.nih.gov/33539890/) published on PubMed.

## Repository Contents

The purpose of the code is to perform a comprehensive network analysis to evaluate the predictive capabilities of various biological networks (such as Funcoup, STRING, and Humannet). It involves constructing and randomizing these networks, selecting seed genes from KEGG pathway data, calculating affinity matrices using a random walk with restart (RWR) algorithm, and assessing performance using precision-recall (PR) and receiver operating characteristic (ROC) curves. The analysis is conducted across multiple iterations and data splits to ensure robustness.
