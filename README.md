# SparseTensorComputations

This repository contains code from my master thesis "Rethinking Sparse Tensor Storage: Incremental Formats as a Path Towards Maximizing Tensor-Vector Multiplication Efficiency" 
It contains the sparse tensor storage format Bit-IF, which is a compression of the coordinate storage format:

Bit-IF: Incremental Sparse Fibres with Bit Encoding

Each mode of the tensor has an increment array, which is the difference between the current and next index. Additionally to the increment, there is a single additional array $b$ of size $d \times nnz_{\mathcal{A}}$ that stores zeros and ones.

