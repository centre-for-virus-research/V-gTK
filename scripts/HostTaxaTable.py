import os
import csv
import sys
import argparse
from os.path import join


class HostTaxaTable:
  def __init__(self,
               gb_matrix,
               output_dir,
               names_file,
               nodes_file,
               base_dir,
               host_output_file,
               child_output_file,
               lineage_output_file):
    self.gb_matrix = gb_matrix
    self.output_dir = output_dir
    self.names_file = names_file
    self.nodes_file = nodes_file
    self.base_dir = base_dir
    self.host_output_file = host_output_file
    self.child_output_file = child_output_file
    self.lineage_output_file = lineage_output_file

  def load_taxa_ids_from_tsv(self):
    column = "host_taxa_id"
    unique = True
    tsv_file = self.gb_matrix
    taxa_ids = []
    with open(tsv_file, newline="", encoding="utf-8") as f:
      reader = csv.DictReader(f, delimiter="\t")
      if column not in reader.fieldnames:
        raise ValueError(f"Column '{column}' not found in TSV file")
      for row in reader:
        raw_val = row.get(column)
        if raw_val is None:
          continue
        val = raw_val.strip()
        if val and val.upper() != "NA":
          taxa_ids.append(int(val))
    if unique:
      return sorted(set(taxa_ids))
    return taxa_ids

  def load_names(self):
    all_names = {}
    sci_names = {}
    with open(self.names_file, encoding="utf-8") as f:
      for line in f:
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 4:
          continue
        taxid_str, name_txt, _, class_name = parts[:4]
        try:
          taxid = int(taxid_str)
        except ValueError:
          continue

        all_names.setdefault(taxid, []).append((name_txt, class_name))
        if class_name == "scientific name":
          sci_names[taxid] = name_txt
    return all_names, sci_names

  def load_nodes(self):
    children_map = {}
    parent_map = {}
    rank_map = {}

    with open(self.nodes_file, encoding="utf-8") as f:
      for line in f:
        parts = [p.strip() for p in line.split("|")]
        if len(parts) < 3:
          continue
        try:
          taxid = int(parts[0])
          parent_taxid = int(parts[1])
        except ValueError:
          continue
        rank = parts[2]

        children_map.setdefault(parent_taxid, []).append(taxid)
        parent_map[taxid] = parent_taxid
        rank_map[taxid] = rank

    return children_map, parent_map, rank_map

  @staticmethod
  def _map_name_type(class_name: str) -> str:
    if class_name == "scientific name":
      return "scientific"
    if class_name in ("genbank common name", "common name"):
      return "generic"
    return "other"

  def _build_lineage_labels(self, taxid, parent_map, rank_map, sci_names):
    lineage = {}
    current = taxid
    visited = set()

    while True:
      if current in visited:
        break
      visited.add(current)

      rank = rank_map.get(current)
      if rank:
        lineage.setdefault(rank, sci_names.get(current, f"taxid_{current}"))

      parent = parent_map.get(current)
      if parent is None or parent == current:
        break
      current = parent

    return lineage

  def write_tables(self):
    os.makedirs(join(self.base_dir, self.output_dir), exist_ok=True)

    children_map, parent_map, rank_map = self.load_nodes()
    all_names, sci_names = self.load_names()
    taxa_list = self.load_taxa_ids_from_tsv()

    host_op_file = join(self.base_dir, self.output_dir, self.host_output_file)
    child_op_file = join(self.base_dir, self.output_dir, self.child_output_file)
    lineage_op_file = join(self.base_dir, self.output_dir, self.lineage_output_file)

    with open(host_op_file, "w", encoding="utf-8", newline="") as out_host:
      writer = csv.writer(out_host, delimiter="\t")
      writer.writerow(["taxa_id", "name", "name_type", "taxonomy_level"])

      for taxid in taxa_list:
        names_for_taxid = all_names.get(taxid, [])
        taxonomy_level = rank_map.get(taxid, "unknown")
        if not names_for_taxid:
          writer.writerow([taxid, "Unknown", "other", taxonomy_level])
        else:
          for name_txt, class_name in names_for_taxid:
            name_type = self._map_name_type(class_name)
            writer.writerow([taxid, name_txt, name_type, taxonomy_level])


    with open(child_op_file, "w", encoding="utf-8", newline="") as out_child:
      writer = csv.writer(out_child, delimiter="\t")
      writer.writerow(["name", "child_taxa_id", "parent_taxa_id",
                       "child_rank", "parent_rank"])

      for parent_taxid in taxa_list:
        children = children_map.get(parent_taxid, [])
        parent_rank = rank_map.get(parent_taxid, "unknown")
        for child_taxid in children:
          child_name = sci_names.get(child_taxid, f"taxid_{child_taxid}")
          child_rank = rank_map.get(child_taxid, "unknown")
          writer.writerow([
            child_name,
            child_taxid,
            parent_taxid,
            child_rank,
            parent_rank
          ])

  
    ranks_order = ["superkingdom", "phylum", "class",
                   "order", "family", "genus", "species"]

    with open(lineage_op_file, "w", encoding="utf-8", newline="") as out_lin:
      writer = csv.writer(out_lin, delimiter="\t")
      writer.writerow(["taxa_id"] + ranks_order)

      for taxid in taxa_list:
        lineage = self._build_lineage_labels(taxid, parent_map, rank_map, sci_names)
        row = [taxid] + [lineage.get(r, "") for r in ranks_order]
        writer.writerow(row)

    print(f"Host taxa table written to {host_op_file}")
    print(f"Child taxa table written to {child_op_file}")
    print(f"Lineage table written to {lineage_op_file}")


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description="Generate host taxa, child taxa and lineage tables from NCBI taxonomy"
  )
  parser.add_argument("-g", "--gb_matrix",
                      default="tmp/GenBank-matrix/gB_matrix_raw.tsv",
                      help="GenBank matrix file")
  parser.add_argument("-o", "--output_dir",
                      default="HostTaxa",
                      help="Output directory (relative to base_dir)")
  parser.add_argument("-n", "--names",
                      default="tmp/Taxa/names.dmp",
                      help="Path to names.dmp file")
  parser.add_argument("-s", "--nodes",
                      default="tmp/Taxa/nodes.dmp",
                      help="Path to nodes.dmp file")
  parser.add_argument("-b", "--base_dir",
                      default="tmp",
                      help="Path to base directory")
  parser.add_argument("-f", "--output_file",
                      default="Host_taxa.tsv",
                      help="Output TSV file for host taxa table")
  parser.add_argument("-c", "--child_output_file",
                      default="Host_taxa_children.tsv",
                      help="Output TSV file for child taxa table")
  parser.add_argument("-y", "--lineage_output_file",
                      default="Host_taxa_lineage.tsv",
                      help="Output TSV file for lineage / hierarchy table")
  args = parser.parse_args()

  write_taxa = HostTaxaTable(
    args.gb_matrix,
    args.output_dir,
    args.names,
    args.nodes,
    args.base_dir,
    args.output_file,
    args.child_output_file,
    args.lineage_output_file,
  )
  write_taxa.write_tables()

