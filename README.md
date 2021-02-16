# VCFmanager

## GWASなど下流の解析に向けたVCFのハンドリングツール

### <解析の流れ>

GATK (<https://github.com/broadinstitute/gatk>) などのバリアントコーラーでVCFを作成  
↓  
00_before_imputation.pyで前処理  
↓  
Beagle (<https://faculty.washington.edu/browning/beagle/beagle.html>)
などのImputationツールで穴埋め  
↓  
10_after_imputation.pyでジェノタイプデータを数値化  
↓  
必要があれば15_transport_txt.shで転置  
↓
