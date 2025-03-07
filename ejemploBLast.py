# blast_example.py - Script básico para realizar búsqueda BLAST usando Biopython
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import os

def blast_sequence(sequence, program="blastn", database="nt", e_value=0.001):
    """
    Realiza una búsqueda BLAST de una secuencia contra la base de datos NCBI
    
    Args:
        sequence: Secuencia en formato string
        program: Tipo de BLAST (blastn, blastp, blastx, etc)
        database: Base de datos para buscar
        e_value: Valor E para filtrar resultados
    
    Returns:
        Los resultados del BLAST en formato XML
    """
    print(f"Ejecutando {program} contra la base de datos {database}...")
    result_handle = NCBIWWW.qblast(program, database, sequence, expect=e_value)
    
    # Guardar resultados en archivo
    with open("blast_results.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    
    result_handle.close()
    print("Búsqueda completada. Resultados guardados en 'blast_results.xml'")
    
    return "blast_results.xml"

def parse_blast_results(xml_file, max_hits=5):
    """
    Analiza los resultados de BLAST y muestra los mejores hits
    
    Args:
        xml_file: Archivo XML con resultados de BLAST
        max_hits: Número máximo de hits a mostrar
    """
    print(f"\nMostrando los {max_hits} mejores resultados:")
    print("-" * 80)
    
    result_handle = open(xml_file)
    blast_records = NCBIXML.parse(result_handle)
    
    for blast_record in blast_records:
        for i, alignment in enumerate(blast_record.alignments):
            if i >= max_hits:
                break
                
            print(f"Secuencia: {alignment.title}")
            for hsp in alignment.hsps:
                print(f"  E-value: {hsp.expect}")
                print(f"  Identidad: {hsp.identities}/{hsp.align_length} ({hsp.identities/hsp.align_length*100:.1f}%)")
                print(f"  Alineamiento: {hsp.query[0:30]}... vs {hsp.sbjct[0:30]}...")
                print()
    
    result_handle.close()

# Secuencia de ejemplo (gen 16S de E. coli)
sequence = """
AGAGTTTGATCCTGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAG
CAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATA
ACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGCACAAAGAGGGGGACCTTAGGGCCTCTTGCCA
TCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTC
TGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATAT
TGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTT
"""

# Ejemplo de uso (comentado para evitar ejecución accidental durante la demostración)
# result_file = blast_sequence(sequence, program="blastn", database="nt", e_value=0.001)
# parse_blast_results(result_file, max_hits=3)

# Nota: Para usar este script, descomenta las últimas dos líneas.
# La primera vez puede tardar dependiendo de la conexión y la carga de NCBI.
