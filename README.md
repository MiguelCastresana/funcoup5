# FunCoup5 Network Analysis

**FunCoup** is a framework to infer genome-wide functional couplings in 21 model organisms. Functional coupling, or functional association, is an unspecific form of association that encompasses direct physical interaction but also more general types of direct or indirect interaction like regulatory interaction or participation in the same process or pathway.

This repository contains the script used to test functional association network performance.

## Background

-  Access the database here: [FunCoup5 website](https://funcoup.org/search/).
-  For a detailed explanation of **FunCoup5** and its applications, please refer to the [FunCoup5 paper](https://pubmed.ncbi.nlm.nih.gov/33539890/).

## Overview

The purpose of the code is to perform a comprehensive network analysis to evaluate the predictive capabilities of various biological networks (such as Funcoup, STRING, and Humannet). It involves constructing and randomizing these networks, selecting seed genes from Orphanet genesets, calculating affinity matrices using a random walk with restart (RWR) algorithm, and assessing performance using precision-recall (PR) and receiver operating characteristic (ROC) curves. The analysis is conducted across multiple iterations and data splits to ensure robustness.

- Loading and filtering gene interaction datasets.
- Processing input gene lists.
- Splitting data for personalized analysis.
- Performing Random Walk with Restart (RWR) on the gene interaction network.
- Calculating and saving performance metrics (ROC and PR curves).

**Contact**:  
Miguel Castresana Aguirre ([miguel.castresana.aguirre@ki.se](mailto:miguel.castresana.aguirre@ki.se))