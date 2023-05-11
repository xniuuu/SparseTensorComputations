# SparseTensorComputations

This repository contains code from my master thesis "Rethinking Sparse Tensor Storage: Incremental Formats as a Path Towards Maximizing Tensor-Vector Multiplication Efficiency" 
It contains the sparse tensor storage format Bit-IF, which is a compression of the coordinate storage format:

**Bit-IF: Incremental Sparse Fibres with Bit Encoding**__
The non-zero values of the tensor are stored in a separate array, with each mode possessing an array that stores the increments which is the difference between the current and next index. This part is inspired from the **I**ncremental **C**ompressed **R**ow **S**torage (ICRS) by Koster. 
Additionally to the increment, there is a single additional array $b$ of size $d \times nnz_{\mathcal{A}}$ that stores zeros and ones. A zero implies no changes in the current mode, while a one prompts an increment. This array $b$ is called the bit encoding array, where each non-zero value is associated with a sequence of bits of size $d$, where $d$ is the dimension size. 

