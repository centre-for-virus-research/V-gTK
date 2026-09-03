"""Repository-level invariants that no single script's tests would catch.

These exist because of two real breakages, both of which passed every unit
test and failed in CI:

1. `scripts/CollectFilteredSequences.py` gained `import projectability` two
   commits before `scripts/projectability.py` was added. Every commit in
   between raised ModuleNotFoundError the moment Nextflow invoked it. Nothing
   imported that script in the test suite, so nothing noticed.

2. The conda environment listed the `defaults` channel. On a runner using the
   preinstalled Anaconda Miniconda that now fails a non-interactive
   Terms-of-Service check, and the whole solve exits 1 with no useful message.

Both are cheap to assert statically and expensive to debug from a CI log.
"""

import ast
import importlib.util
import re
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = REPO_ROOT / "scripts"
WORKFLOWS_DIR = REPO_ROOT / ".github" / "workflows"

SCRIPT_FILES = sorted(SCRIPTS_DIR.glob("*.py"))

#: Scripts known to import something that is not in the repository, with the
#: reason. Entries are asserted to STILL be broken, so fixing one fails this
#: test and tells you to delete the entry - a known gap stays visible instead
#: of quietly becoming permanent.
#:
#: Nothing here may be reachable from the pipeline; a broken script that
#: vgtk-init.nf invokes is a build failure, not a known gap.
KNOWN_UNRESOLVED = {
	"HcvBackgroundNotebook.py": (
		"imports PlotHcvMutationTrends, which has never been committed. The "
		"script also hardcodes /home3/... paths, so it appears to be a personal "
		"analysis script from the other development machine. Not referenced by "
		"vgtk-init.nf, vgtk-rabv.sh or nextflow.config."
	),
}


def _declared_channels(text: str):
    """Pull the `channels:` list out of environment.yml.

    Parsed by hand rather than with PyYAML: this file guards the environment,
    so it must not itself depend on something the environment might not have.
    """
    channels, inside = [], False
    for line in text.splitlines():
        if re.match(r"^channels:\s*$", line):
            inside = True
            continue
        if inside:
            if re.match(r"^\s*#", line):
                continue
            item = re.match(r"^\s*-\s*(\S+)", line)
            if item:
                channels.append(item.group(1))
                continue
            if line.strip():          # a new top-level key ends the list
                break
    return channels


def _imported_top_level_names(path: Path):
    """Top-level module names a file imports, ignoring relative imports."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                names.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.level:          # relative import, not our concern
                continue
            if node.module:
                names.add(node.module.split(".")[0])
    return names


def _is_resolvable(name: str) -> bool:
    """A name resolves if it is a sibling script or is installed."""
    if (SCRIPTS_DIR / f"{name}.py").exists():
        return True
    if (SCRIPTS_DIR / name / "__init__.py").exists():
        return True
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError, ModuleNotFoundError):
        return False


class TestEveryImportResolves:
    """The projectability guard.

    A script that imports a module which does not exist is broken for every
    caller, but only fails when something actually runs it - which, for a
    pipeline script, means in the middle of a Nextflow task.
    """

    def test_there_are_scripts_to_check(self):
        # A glob that silently matches nothing would make this file vacuous.
        assert len(SCRIPT_FILES) > 30

    @pytest.mark.parametrize("script", SCRIPT_FILES, ids=lambda p: p.name)
    def test_every_imported_module_exists(self, script):
        unresolved = sorted(
            name for name in _imported_top_level_names(script)
            if not _is_resolvable(name)
        )
        if script.name in KNOWN_UNRESOLVED:
            assert unresolved, (
                f"{script.name} is listed in KNOWN_UNRESOLVED but all its "
                f"imports now resolve - delete the entry"
            )
            return
        assert not unresolved, (
            f"{script.name} imports {unresolved}, which is neither a sibling "
            f"module in scripts/ nor installed in this environment"
        )

    @pytest.mark.parametrize("name", sorted(KNOWN_UNRESOLVED))
    def test_known_gaps_are_not_reachable_from_the_pipeline(self, name):
        """A broken script the pipeline calls is a failure, not a known gap."""
        stem = Path(name).stem
        for caller in ("vgtk-init.nf", "vgtk-rabv.sh", "nextflow.config"):
            path = REPO_ROOT / caller
            if path.exists():
                assert stem not in path.read_text(encoding="utf-8"), (
                    f"{name} is exempted from the import check but {caller} "
                    f"invokes it - fix the script instead of exempting it"
                )


class TestCondaChannels:
    """The `defaults` guard.

    Reinstating it breaks the Nextflow CI job and reports only
    "conda failed with exit code 1", which is a long way from the cause.
    """

    def test_environment_yml_does_not_use_defaults(self):
        channels = _declared_channels((REPO_ROOT / "environment.yml").read_text())
        assert "defaults" not in channels, (
            "the `defaults` channel fails conda's non-interactive "
            "Terms-of-Service check on CI runners"
        )

    def test_environment_yml_keeps_the_channels_it_needs(self):
        channels = _declared_channels((REPO_ROOT / "environment.yml").read_text())
        assert "conda-forge" in channels and "bioconda" in channels

    def test_the_channel_parser_actually_finds_them(self):
        # Guards the guard: a parser that silently returned [] would make the
        # assertion above pass for any file at all.
        assert _declared_channels(
            "name: x\nchannels:\n  # a comment\n  - conda-forge\n  - bioconda\n"
            "dependencies:\n  - python\n"
        ) == ["conda-forge", "bioconda"]

    @pytest.mark.parametrize(
        "workflow", sorted(WORKFLOWS_DIR.glob("*.yml")), ids=lambda p: p.name
    )
    def test_no_workflow_asks_for_defaults(self, workflow):
        for number, line in enumerate(workflow.read_text().splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue                      # commented-out jobs are inert
            if stripped.startswith("channels:"):
                assert "defaults" not in stripped, f"{workflow.name}:{number}"


class TestWorkflowsAreParseable:
    @pytest.mark.parametrize(
        "workflow", sorted(WORKFLOWS_DIR.glob("*.yml")), ids=lambda p: p.name
    )
    def test_workflow_yaml_is_valid(self, workflow):
        # A workflow that does not parse simply never runs, and GitHub reports
        # that as the job not existing rather than as an error. PyYAML is not a
        # pipeline dependency, so this checks only when it happens to be there.
        yaml = pytest.importorskip("yaml")
        assert yaml.safe_load(workflow.read_text()) is not None
