"""Submit prediction results to CASP by email
"""
import os
from email.message import EmailMessage
import gzip
from io import StringIO
import smtplib

from Bio.PDB import PDBParser, MMCIFParser, PDBIO

def yield_chains(input_file):
  input_pid = os.path.basename(input_file)
  if input_file.endswith(".gz"):
    input_pid, _ = os.path.splitext(input_pid)
  input_pid, input_type = os.path.splitext(input_pid)
  assert input_type in (".cif", ".pdb")
  pdb_parser = PDBParser(
      QUIET=True) if input_type in (".pdb",) else MMCIFParser(QUIET=True)

  if input_file.endswith(".gz"):
    with gzip.open(input_file, "rt") as f:
      protein_structure = pdb_parser.get_structure("none", f)
  else:
    protein_structure = pdb_parser.get_structure("none", input_file)
  model = list(protein_structure.get_models())[0]
  io = PDBIO()
  io.set_structure(protein_structure)
  for chain in model.get_chains():
    f = StringIO()
    io.save(f)
    yield chain.id, f.getvalue()

def main(args):  # pylint: disable=redefined-outer-name
  mail_text = (f"PFRMAT {args.format}\n"
               f"TARGET {args.target}\n"
               f"AUTHOR {args.author}\n"
                "REMARK FALCON2\n"
               f"METHOD {args.method}\n"
               f"MODEL  {args.model}\n"
               f"STOICH {args.stoich}\n")

  for input_file in args.input_file:
    for _, input_text in yield_chains(input_file):
      if args.format == "TS":
        input_text = filter(lambda x: x.startswith("ATOM"),
                            input_text.splitlines())
        input_text = "\n".join(input_text)

      mail_text += f"PARENT {args.parent}\n"
      mail_text += f"{input_text}\n"
      mail_text +=  "TER\n"

  mail_text +=  "END"

  msg = EmailMessage()
  msg["Subject"] = f"{args.target} {args.format} {args.method}@{args.model}"
  msg["From"] = args.sender
  msg["To"] = args.receiver
  msg.add_header("Content-Type", "text/plain")
  # msg.add_header("Content-Transfer-Encoding", "base64")
  msg.set_payload(mail_text)
  if args.verbose:
    print(msg.as_string())

  if args.output_file:
    with open(args.output_file, "w") as f:
      f.write(mail_text)
  if not args.dry_run:
    s = smtplib.SMTP_SSL("mail.cstnet.cn")
    if args.passwd is not None:
      s.login(args.sender, args.passwd)  # "ZBZBHYFG8DC*Bjg!")
    r = s.send_message(msg)
    s.quit()
  if args.verbose:
    print("---------------")
    print(mail_text)
    print("---------------")
  if not args.dry_run:
    print(f"Done: {r}")

if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument("input_file", type=str, nargs="+", help="pdb/cif file.")
  parser.add_argument("--output_file", type=str, default=None,
                      help="output file.")
  parser.add_argument("--sender", type=str, default="falcon@ict.ac.cn",
                      help="email sender.")
  parser.add_argument("--passwd", type=str, default=None,
                      help="sender password.")
  parser.add_argument("--receiver", type=str,
                      default="models@predictioncenter.org",
                      help="email receiver.")
  parser.add_argument("--format", type=str, default="TS",
                      choices=["TS", "QA", "LG"],
                      help="format specification code: TS, QA, LG")
  parser.add_argument("--target", type=str, required=True,
                      help="target name.")
  parser.add_argument("--author", type=str, default="5425-8237-3967",
                      help="registration code of server group name.")
  parser.add_argument("--method", type=str, default="FALCON2",
                      help="records describing the methods used.")
  parser.add_argument("--stoich", type=str, default="A2",
                      help="stoichiometry information for each target.")
  parser.add_argument("--model", type=int, default=1,
                      help="model index.")
  parser.add_argument("--parent", type=str, default="N/A",
                      help="parent.")
  parser.add_argument("--dry_run", action="store_true",
                      help="dry run.")
  parser.add_argument("-v", "--verbose", action="store_true",
                      help="verbose.")

  args = parser.parse_args()

  main(args)
