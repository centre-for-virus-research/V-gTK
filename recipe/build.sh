#!/usr/bin/env bash

set -euo pipefail

# Create V-gTK installation directories
mkdir -p "${PREFIX}/share/vgtk/scripts"
mkdir -p "${PREFIX}/bin"

# Copy all existing V-gTK scripts
cp -R "${SRC_DIR}/scripts/"* \
    "${PREFIX}/share/vgtk/scripts/"

# Copy documentation
if [ -f "${SRC_DIR}/README.md" ]; then
    cp "${SRC_DIR}/README.md" "${PREFIX}/share/vgtk/"
fi

if [ -f "${SRC_DIR}/LICENSE" ]; then
    cp "${SRC_DIR}/LICENSE" "${PREFIX}/share/vgtk/"
fi

# Create the V-gTK launcher
cat > "${PREFIX}/bin/vgtk" <<'EOF'
#!/usr/bin/env python3

import os
import sys
import subprocess


def main():

    if len(sys.argv) < 2:
        print("V-gTK")
        print("")
        print("Usage:")
        print("  vgtk <command> [arguments]")
        print("")
        print("Examples:")
        print("  vgtk DownloadGFF --help")
        print("  vgtk GenBankParser --help")
        sys.exit(0)

    command = sys.argv[1]

    if command.endswith(".py"):
        command = command[:-3]

    scripts_dir = os.path.join(
        sys.prefix,
        "share",
        "vgtk",
        "scripts"
    )

    script_path = os.path.join(
        scripts_dir,
        command + ".py"
    )

    if not os.path.isfile(script_path):
        print(f"V-gTK command not found: {command}")
        sys.exit(1)

    cmd = [
        sys.executable,
        script_path
    ] + sys.argv[2:]

    result = subprocess.run(cmd)

    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
EOF

chmod +x "${PREFIX}/bin/vgtk"
