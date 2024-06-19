
def load_datasets():
    """Load dataset locations. Update to filelist?"""
    datasets={'QM9':"qm9_data.dump",
              'Quinones':"tabor_nosulf.dump",
              'Wang':'wang_data.dump',
              'Gecko':'gecko_full.dump',
              'ExpMoNA':"allmona_cleaned.dump",
              'nablaDFT':"nablaDFT.dump",
             'mb_eu':"mb_eu_cleaned.dump"}
    return datasets
