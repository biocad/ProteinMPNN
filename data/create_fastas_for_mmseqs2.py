from pathlib import Path
from enum import Enum
from tqdm import tqdm
import logging
import sys
import traceback
from pathlib import Path
from typing import Iterable, cast
from argparse import ArgumentParser
import pandas as pd

from proteinlib.structure.antibody_antigen_complex import AntibodyAntigenComplex, NumberingScheme

class Regions(Enum):
    FR_L1=0
    CDR_L1=1
    FR_L2=2
    CDR_L2=3
    FR_L3=4
    CDR_L3=5
    FR_L4=6

    FR_H1=7
    CDR_H1=8
    FR_H2=9
    CDR_H2=10
    FR_H3=11
    CDR_H3=12
    FR_H4=13

if __name__=='__main__':
    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    
    parser=ArgumentParser('Create fastas for clustering via MMseqs')
    parser.add_argument('--summary_csv', type=str,default=Path('/mnt/sabdab/summary.csv'))
    parser.add_argument('--chothia_subdir', type=str,default=Path("/mnt/sabdab/chothia"))
    parser.add_argument('--regions',default=['CDR_H3','CDR_H2','CDR_H1','CDR_L1','CDR_L2','CDR_L3'],nargs='*')

    args=parser.parse_args()

    summary_csv=args.summary_csv
    chothia_subdir=args.chothia_subdir
    regions_list=args.regions

    # only protein and peptide antigens
    df = (
        pd.read_csv(
            summary_csv,
            sep="\t",
            usecols=["pdb", "Hchain", "Lchain", "antigen_chain", "antigen_type"],
        )
        .query("antigen_type in ('protein', 'peptide')")
        .dropna()
        .drop_duplicates('pdb')
        .reset_index()
    )
    print(f"Summary records: {df.shape[0]}")

    d={}
    for i, row in tqdm(cast(Iterable[tuple[int, pd.Series]], df.iterrows()),total=df.shape[0]):
        uid = f'{row["pdb"]}_{row["Hchain"]}+{row["Lchain"]}-{row["antigen_chain"]}'
        try:
            sequences: list[tuple[str, str]] = []
            # row = cast(pd.Series, row)
            antigen_chains = tuple(map(lambda s: s.strip(), str(row["antigen_chain"]).split(" | ")))
            ab_complex = AntibodyAntigenComplex.from_pdb(
                pdb=chothia_subdir / f"{row['pdb']}.pdb",
                heavy_chain_id=str(row["Hchain"]),
                light_chain_id=str(row["Lchain"]),
                antigen_chain_ids=antigen_chains,
                numbering=NumberingScheme.CHOTHIA,
            )
            # get regions
            vh_regions = [region.sequence for region in list(ab_complex.antibody.heavy_chain.regions)]
            vl_regions = [region.sequence for region in list(ab_complex.antibody.light_chain.regions)]
            all_regions = [*vl_regions, *vh_regions]

            for region in regions_list:
                sequences.append(all_regions[getattr(Regions,region).value])
            d[uid]=''.join(sequences)
        except FileNotFoundError:
            continue
        except Exception as err:
            logging.warning(f"In complex {uid}: {traceback.format_exception(err)}")
            continue
    fasta_str=''
    for k,v in d.items():
        fasta_str='\n'.join([fasta_str,f">{k}\n{v}"])

    fasta_str=fasta_str.strip()
    output_path=Path(f'DB_{"_".join(regions_list)}.fasta')
    output_path.write_text(fasta_str)

# after creating fast launch clustering via MMseqs2
# FASTA=<fasta file name> && IDENTITY=<identity parameter for clustering> && mmseqs easy-cluster $FASTA  clusterRes_${IDENTITY}_${FASTA}  tmp --min-seq-id $IDENTITY -c 0.8 --cov-mode 1
# example
# FASTA=DB_FR_H1_FR_H2_FR_H3_FR_L1_FR_L2_FR_L3.fasta && IDENTITY=0.5 && mmseqs easy-cluster $FASTA  clusterRes_${IDENTITY}_${FASTA}  tmp --min-seq-id $IDENTITY -c 0.8 --cov-mode 1