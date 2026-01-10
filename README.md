# Chou-Fasman Algorithm for Protein Secondary Structure Prediction

A Python implementation of the Chou-Fasman algorithm for predicting protein secondary structure from amino acid sequences.

## 🧬 Overview

The Chou-Fasman algorithm is a classic method for predicting protein secondary structure (α-helices and β-sheets) based on the propensity of amino acids to form these structures. This implementation provides accurate predictions with conflict resolution and detailed analysis.

## ✨ Features

- **Secondary Structure Prediction**: Predicts α-helices and β-sheets from protein sequences
- **Propensity-Based Analysis**: Uses empirically derived propensity values for each amino acid
- **Conflict Resolution**: Handles overlapping predictions intelligently
- **Detailed Output**: Provides comprehensive analysis including:
  - Initial predictions with conflicts
  - Final predictions after conflict resolution
  - Sequence alignment visualization
  - Statistical summaries

## 🔬 Algorithm Details

### Prediction Rules

**α-Helix Prediction:**
- Window size: 6 amino acids
- Threshold: ≥4 residues with P(α) ≥ 1.0
- Extension: Continue while average P(α) ≥ 1.0 in 4-residue windows

**β-Sheet Prediction:**
- Window size: 5 amino acids  
- Threshold: ≥3 residues with P(β) ≥ 1.0
- Extension: Continue while average P(β) ≥ 1.0 in 4-residue windows

**Conflict Resolution:**
- Compare sum of propensities for conflicting regions
- Assign structure type with higher total propensity

## 📊 Propensity Values

The algorithm uses experimentally determined propensity values:

**High α-helix formers:** Glu (1.53), Ala (1.45), Leu (1.34)
**High β-sheet formers:** Met (1.67), Val (1.65), Ile (1.60)
**Structure breakers:** Pro (α: 0.59, β: 0.62), Gly (α: 0.53, β: 0.81)

## 🚀 Usage

### Basic Usage
```python
python chou_fasman_algorithm.py
```

### Input Format
The algorithm accepts single-letter amino acid sequences:
```python
sequence = "MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQDWAYAASKESHATLVFHNLLGEIDQQYSRFLQESNVLYQHNLRRIKQFLQSRYLEKPMEIARIVARCLWEESRLLQTAATAAQQGGQANHPTAAVVTEKQQMLEQHLQDVRKRVQDLEQKMKVVENLQDDFDFNYKTLKSQGDMQDLNGNNQSVTRQKMQQLEQMLTALDQMRRSIVSELAGLLSAMEYVQKTLTDEELADWKRRQQIACIGGPPNICLDRLENWITSLAESQLQTRQQIKKLEELQQKVSYKGDPIVQHRPMLEERIVELFRNLMKSAFVVERQPCMPMHPDRPLVIKTGVQFTTKVRLLVKFPELNYQLKIKVCIDKDSGDVAALRGSRKFNILGTNTKVMNMEESNNGSLSAEFKHLTLREQRCGNGGRANCDASLIVTEELHLITFETEVYHQGLKIDLETHSLPVVVISNICQMPNAWASILWYNMLTNNPKNVNFFTKPPIGTWDQVAEVLSWQFSSTTKRGLSIEQLTTLAEKLLGPGVNYSGCQITWAKFCKENMAGKGFSFWVWLDNIIDLVKKYILALWNEGYIMGFISKERERAILSTKPPGTFLLRFSESSKEGGVTFTWVEKDISGKTQIQSVEPYTKQQLNNMSFAEIIMGYKIMDATNILVSPLVYLYPDIPKEEAFGKYCRPESQEHPEADPGSAAPYLKTKFICVTPTTCSNTIDLPMSPRTLDSLMQFGNNGEGAEPSAGGQFESLTFDMELTSECATSPM"
```

## 📈 Output Format

The program provides detailed output including:

1. **Initial Predictions** (with conflicts)
2. **Conflict Analysis** with resolution details
3. **Final Structure Assignment**
4. **Sequence Alignment** showing predicted structures

Example output:
```
Helix regions (including conflicted regions):
Helix: positions (15-28) -> YLEQLHQLYSDSFPM

Beta strand regions (including conflicted regions):
Beta strand: positions (25-35) -> SFPMELRQFLA

Conflicting region:
Conflict positions(25-28) -> seq = SFPM
sum P_alpha = 4.12, sum P_beta = 5.23 => assigned to S

Final full-sequence alignment:
Sequence : MAQWNQLQQLDTRYLEQLHQLYSDSFPMELRQFLAPWIESQ...
Structure: ----HHHHHHHHHHHHHH----SSSSSSSSSS---------...
```

## 🧪 Algorithm Validation

This implementation follows the original Chou-Fasman methodology:
- Uses published propensity values
- Implements standard window sizes and thresholds
- Includes proper extension and conflict resolution rules
- Provides 0-based indexing for computational accuracy

## 📋 Requirements

- Python 3.6+
- No external dependencies required

## 🔧 Technical Details

- **Input Processing**: Converts single-letter to three-letter amino acid codes
- **Window Scanning**: Systematic analysis of sequence windows
- **Extension Algorithm**: Bidirectional extension based on average propensities
- **Conflict Resolution**: Quantitative comparison of competing predictions
- **Output Formatting**: Clear visualization of results

## 📚 Scientific Background

The Chou-Fasman algorithm was developed by Peter Y. Chou and Gerald D. Fasman in 1974. It was one of the first computational methods for protein secondary structure prediction and remains valuable for:

- Educational purposes in bioinformatics
- Rapid preliminary structure analysis
- Comparison with modern prediction methods
- Understanding structure-sequence relationships

## 🤝 Contributing

Contributions are welcome! Areas for improvement:
- Additional validation datasets
- Performance optimizations
- Extended output formats
- Integration with modern prediction methods

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact

**Author**: Arnav Garg  
**Student ID**: 2024107  
**Institution**: IIIT Delhi

## 🙏 Acknowledgments

- Original algorithm by Chou & Fasman (1974)
- Propensity values from empirical studies
- Bioinformatics community for continued development

---

*This implementation is part of computational biology coursework and demonstrates classical approaches to protein structure prediction.*