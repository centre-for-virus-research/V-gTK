import os
import csv
import time
import subprocess
from os.path import join
from argparse import ArgumentParser


class SoftwareVersionChecker:
    def __init__(self, tmp_dir, output_dir, table_name, environment_yml):
        self.tmp_dir = tmp_dir
        self.output_dir = output_dir
        self.table_name = table_name
        self.environment_yml = environment_yml

        self.software_commands = {
            "Nextalign": ["nextalign", "--version"],
            "BLAST": ["blastn", "-version"],
            "MAFFT": ["mafft", "--version"],
            "IQ-TREE": ["iqtree2", "--version"],
            "FastTree": ["FastTree", "-help"],
            "Python": ["python", "--version"]
        }

        self.software_versions = {}
        self.environment_name = ""
        self.environment_channels = []
        self.environment_dependencies = []

        self.check_software_versions()
        self.read_environment_yml()
        self.save_versions_to_tsv()

    def get_software_version(self, software, command):
        try:
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            output = result.stdout if result.stdout else result.stderr

            lines = output.split("\n")
            for line in lines:
                if "version" in line.lower() or software.lower() in line.lower():
                    return line.strip()

            return output.strip().split("\n")[0]

        except FileNotFoundError:
            return "Not Installed"

    def check_software_versions(self):
        self.software_versions = {
            software: self.get_software_version(software, command)
            for software, command in self.software_commands.items()
        }

    def split_dependency_version(self, dependency):
        dependency = dependency.strip().strip('"').strip("'")

        if "==" in dependency:
            package, version = dependency.split("==", 1)
            return package.strip(), version.strip()

        if "=" in dependency:
            package, version = dependency.split("=", 1)
            return package.strip(), version.strip()

        return dependency.strip(), ""

    def read_environment_yml(self):
        if not self.environment_yml:
            return

        if not os.path.exists(self.environment_yml):
            self.environment_dependencies.append({
                "source": "environment",
                "package": "environment.yml",
                "version": f"Not found: {self.environment_yml}"
            })
            return

        current_section = None
        in_pip_section = False

        with open(self.environment_yml, "r") as file:
            for raw_line in file:
                line = raw_line.rstrip("\n")
                stripped = line.strip()

                if not stripped or stripped.startswith("#"):
                    continue

                if stripped.startswith("name:"):
                    self.environment_name = stripped.split(":", 1)[1].strip()
                    current_section = None
                    in_pip_section = False
                    continue

                if stripped == "channels:":
                    current_section = "channels"
                    in_pip_section = False
                    continue

                if stripped == "dependencies:":
                    current_section = "dependencies"
                    in_pip_section = False
                    continue

                if current_section == "channels":
                    if stripped.startswith("- "):
                        channel = stripped[2:].strip()
                        self.environment_channels.append(channel)
                    continue

                if current_section == "dependencies":
                    if stripped.startswith("- pip:"):
                        in_pip_section = True
                        continue

                    if stripped.startswith("- "):
                        dependency = stripped[2:].strip()

                        package, version = self.split_dependency_version(dependency)

                        self.environment_dependencies.append({
                            "source": "conda",
                            "package": package,
                            "version": version
                        })

                        in_pip_section = False
                        continue

                    if in_pip_section and stripped.startswith("- "):
                        dependency = stripped[2:].strip()

                        package, version = self.split_dependency_version(dependency)

                        self.environment_dependencies.append({
                            "source": "pip",
                            "package": package,
                            "version": version
                        })

    def save_versions_to_tsv(self):
        os.makedirs(self.tmp_dir, exist_ok=True)
        os.makedirs(join(self.tmp_dir, self.output_dir), exist_ok=True)

        output_path = os.path.join(self.tmp_dir, self.output_dir, self.table_name)

        with open(output_path, "w", newline="") as file:
            writer = csv.writer(file, delimiter="\t")

            writer.writerow(["Software", "Version"])

            for software, version in self.software_versions.items():
                writer.writerow([software, version])

            writer.writerow(["Time of creation", time.time()])
            writer.writerow(["vgtk version", "v.1.0.0"])

            writer.writerow([])
            writer.writerow(["Environment item", "Version"])

            if self.environment_yml:
                writer.writerow(["Environment file", self.environment_yml])

            if self.environment_name:
                writer.writerow(["Environment name", self.environment_name])

            if self.environment_channels:
                writer.writerow(["Environment channels", ", ".join(self.environment_channels)])

            for dependency in self.environment_dependencies:
                source = dependency["source"]
                package = dependency["package"]
                version = dependency["version"]

                writer.writerow([
                    f"{source}: {package}",
                    version
                ])


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Create a TSV with available software version information"
    )

    parser.add_argument(
        "-d",
        "--tmp_dir",
        help="Temp directory",
        default="tmp"
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        help="Output directory",
        default="Software_info"
    )

    parser.add_argument(
        "-f",
        "--table_name",
        help="Name of the TSV to be saved",
        default="software_info.tsv"
    )

    parser.add_argument(
        "-e",
        "--environment_yml",
        help="Path to environment.yml file",
        default="environment.yml"
    )

    args = parser.parse_args()

    checker = SoftwareVersionChecker(
        args.tmp_dir,
        args.output_dir,
        args.table_name,
        args.environment_yml
    )