![CFDPSM](https://github.com/user-attachments/assets/bb1c20b9-2e3a-48c1-bd47-e3624951b4d3)


# 1. EECFS

> In this section, we propose the EECFS algorithm. During the spouse discovery phase, the EDMB algorithm identifies the PC set of each variable in the target's PC as candidate spouses. However, spouse nodes exist only in certain children's PC sets, causing EDMB to perform many unnecessary CI tests for all parents and children, which increases computational overhead. To address this, EECFS directly selects spouse candidates from the PC sets of target children with multiple parents, reducing CI tests and improving efficiency.

# 2. CFDPSM

> CFDPSM consists of three sequential modules:
I. Feature Annotation: Features are extracted at the DNA, RNA, and protein levels to annotate each variant. Features with a missing value ratio ≤ 5% are retained, and missing values are imputed using MissForest\cite{stekhoven2012missforest}.
II. Causal Feature Selection: For the raw feature set, multiple causal feature selection methods are applied separately to identify their respective MB sets. From the extracted MB candidates, the subset with the minimal dimensionality that achieves optimal or near-optimal performance on the training set is selected as the final feature subset.
III. Classification: The MB feature subsets from the three molecular levels are concatenated to form the final feature set. LightGBM is then employed for model training and prediction. sSNVs with predicted scores > 0.5 are classified as pathogenic, while those with scores ≤ 0.5 are classified as benign.

# 3. Directory Structure

> NOTE: Dataset are available in FigShare: 10.6084/m9.figshare.31566328.

- 📂 code/: Main code directory
	- 📂 Causal_feature_selection/
		- 📂 alg_MB/: Other competing causal feature selection methods
			- 📄 \_G2 suffix indicates G² test, \_Z suffix indicates Fisher Z test
		- 📂 Benchmark_NB_dataset/: Benchmark Bayesian network datasets for comparison
		- 📂 common/: General MATLAB function files
		- 📂 evaluation/: MATLAB evaluation scripts
		- 📂 Real_world_dataset/: Real-world datasets
			- 📄 make_data.m: Generates 10-fold cross-validation index data
		- 📂 n_XXXX/: n corresponds to the dataset number, XXXX corresponds to the dataset name (see Supplementary Table S2)
		- 📄 example_MB.m: Demonstrates Markov Blanket learning using the Causal Learner framework
		- 📄 Causal_Learner.m: Implements a unified interface for causal structure and Markov Blanket learning from discrete or continuous data
		- 📄 extract_MB.py: Extracts MB feature data
	- 📄 5fold_lightgbm.py: Performs 5-fold training using the best parameters
	- 📄 call_5fold_lightgbm.py: Calls the 5fold_lightgbm.py script
	- 📄 lightgbm_test.py: Evaluates the final model
- 📂 data/: Training and evaluation data
	- 📂 DNA/: DNA-level feature space for synonymous variant prediction task; train_ and test_ prefixes indicate training and test sets; missing values - imputed via MissForest
	- 📂 RNA/: RNA-level feature space, same structure as DNA
	- 📂 Protein/: Protein-level feature space, same structure as DNA
	- 📂 MB_feature/: Causal feature sets selected at three molecular levels using different causal feature selection methods
- 💻 model/: Final model weights of CFDPSM

# 4. Citation Information

1. Other MB methods in this repository (except EDMB) are provided via the Causal Learner package: http://bigdata.ahu.edu.cn/causal-learner. Please cite the following if used:

```bibtex
@article{ling2022causal,
  title={Causal learner: A toolbox for causal structure and markov blanket learning},
  author={Ling, Zhaolong and Yu, Kui and Zhang, Yiwen and Liu, Lin and Li, Jiuyong},
  journal={Pattern Recognition Letters},
  volume={163},
  pages={92--95},
  year={2022},
  publisher={Elsevier}
}
```

2. EECFS & CFDPSM

```bibtex
In preparation.
```
