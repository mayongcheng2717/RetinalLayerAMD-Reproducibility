from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "Additional_Files/Additional_file_1_Supplementary_Tables_S1-S13.xlsx"
REPLACEMENTS = {
    "• This workbook accompanies the final manuscript and contains complete result tables rather than a reviewer-response document.": "• This workbook accompanies the final manuscript and contains complete result tables; it is not a reviewer-response document.",
    "• Source IDs map to S13a_data_sources and to the reproducibility archive. No public code repository or archived DOI is claimed.": "• Source IDs map to S13a_data_sources and to the reproducibility archive. The GitHub repository location is specified in the Availability of data and materials statement; no archived DOI has been assigned.",
}

def update_workbook(path=WORKBOOK):
    path = Path(path)
    with zipfile.ZipFile(path, "r") as z:
        shared = z.read("xl/sharedStrings.xml").decode("utf-8")
    if all(old not in shared and new in shared for old, new in REPLACEMENTS.items()):
        print("Workbook README is already current; no file was changed.")
        return
    counts = {old: shared.count(old) for old in REPLACEMENTS}
    if counts != {old: 1 for old in REPLACEMENTS}:
        raise RuntimeError(f"Unexpected README text state: {counts}")
    temporary = path.with_suffix(".tmp.xlsx")
    with zipfile.ZipFile(path, "r") as zin, zipfile.ZipFile(temporary, "w") as zout:
        for info in zin.infolist():
            data = zin.read(info.filename)
            if info.filename == "xl/sharedStrings.xml":
                text = data.decode("utf-8")
                for old, new in REPLACEMENTS.items():
                    text = text.replace(old, new)
                data = text.encode("utf-8")
            zout.writestr(info, data)
    temporary.replace(path)
    with zipfile.ZipFile(path) as z:
        assert z.testzip() is None
    print("Updated only xl/sharedStrings.xml; formulas, cached values, and style parts were not edited.")

if __name__ == "__main__":
    update_workbook()
