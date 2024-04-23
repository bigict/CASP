"""Submit prediction results to CASP by email
"""
from email.message import EmailMessage
import smtplib

def main(args):  # pylint: disable=redefined-outer-name
  with open(args.input_file, "r") as f:
    input_text = f.read()

  if args.format == "TS":
    input_text = filter(lambda x: x.startswith("ATOM"), input_text.splitlines())
    input_text = "\n".join(input_text)

  input_text = (f"PFRMAT {args.format}\n"
                f"TARGET {args.target}\n"
                f"AUTHOR {args.author}\n"
                f"REMARK FALCON2\n"
                f"METHOD {args.method}\n"
                f"MODEL  {args.model}\n"
                f"PARENT {args.parent}\n"
                f"{input_text}\n"
                f"TER\n"
                f"END")
  msg = EmailMessage()
  msg["Subject"] = f"{args.target} {args.format} {args.method}@{args.model}"
  msg["From"] = args.sender
  msg["To"] = args.receiver
  msg.add_header("Content-Type", "text/plain")
  # msg.add_header("Content-Transfer-Encoding", "base64")
  msg.set_payload(input_text)
  if args.verbose:
    print(msg.as_string())

  s = smtplib.SMTP_SSL("mail.cstnet.cn")
  if args.passwd is not None:
    s.login(args.sender, args.passwd)  # "ZBZBHYFG8DC*Bjg!")
  r = s.send_message(msg)
  s.quit()
  if args.verbose:
    print("---------------")
    print(input_text)
    print("---------------")
  print(f"Done: {r}")

if __name__ == "__main__":
  import argparse

  parser = argparse.ArgumentParser(
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument("input_file", type=str, help="result file.")
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
  parser.add_argument("--model", type=int, default=1,
                      help="model index.")
  parser.add_argument("--parent", type=str, default="N/A",
                      help="parent.")
  parser.add_argument("-v", "--verbose", action="store_true",
                      help="verbose")

  args = parser.parse_args()

  main(args)
