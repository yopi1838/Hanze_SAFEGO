
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# USER PARAMETERS
# ============================================================

INCLUDE_TOP_STEEL_BEAM = False

G = 9.8104
WALL_HEIGHT = 2.58       # m
WALL_WEIGHT = 12.94      # kN
SLAB1_MASS = 165.0       # kg
SLAB2_MASS = 165.0       # kg
STEEL_BEAM_MASS = 61.4   # kg

TESTS = {
    9: "US1",
    10: "ST1-FDM",
    11: "ST2-HB",
    12: "US2"
}

OUTPUT = Path("BaseShear_Plots")
OUTPUT.mkdir(exist_ok=True)


# ============================================================
# FUNCTIONS
# ============================================================

def read_heights(test_number):
    info_file = f"Test{test_number}/Test{test_number}_Info.xlsx"
    info = pd.read_excel(info_file, header=None)

    # F16 -> acceleration column 15
    # F17 -> acceleration column 16
    # F18 -> acceleration column 17
    h1 = float(info.iloc[15, 5])
    h2 = float(info.iloc[16, 5])
    h3 = float(info.iloc[17, 5])

    if not (0 < h1 < h2 < h3 < WALL_HEIGHT):
        raise ValueError(f"Invalid heights: {h1}, {h2}, {h3}")

    return h1, h2, h3


def calculate_base_shear(df, heights):
    h1, h2, h3 = heights
    H = WALL_HEIGHT

    # Tributary height fractions
    f1 = (h1 + (h2 - h1) / 2) / H
    f2 = ((h2 - h1) / 2 + (h3 - h2) / 2) / H
    f3 = (H - h3 + (h3 - h2) / 2) / H

    # Accelerations in g
    a13 = pd.to_numeric(df["13"], errors="coerce")
    a14 = pd.to_numeric(df["14"], errors="coerce")
    a15 = pd.to_numeric(df["15"], errors="coerce")
    a16 = pd.to_numeric(df["16"], errors="coerce")
    a17 = pd.to_numeric(df["17"], errors="coerce")

    # Slab inertia forces (kN)
    V_slabs = G / 1000 * (
        SLAB1_MASS * a13 +
        SLAB2_MASS * a14
    )

    # Wall inertia forces (kN)
    V_wall = WALL_WEIGHT * (
        f1 * a15 +
        f2 * a16 +
        f3 * a17
    )

    V = V_slabs + V_wall

    # Optional top steel beam inertia
    if INCLUDE_TOP_STEEL_BEAM:
        V += G / 1000 * STEEL_BEAM_MASS * a17

    return V


# ============================================================
# PROCESS ALL TESTS
# ============================================================

for test_number, specimen in TESTS.items():

    folder = Path(f"Test{test_number}")
    heights = read_heights(test_number)

    print(f"{specimen}: accelerometer heights = {heights}")

    files = sorted(
        folder.glob(f"Test{test_number}Run*.xlsx"),
        key=lambda p: int(p.stem.split("Run")[-1])
    )

    fig, ax = plt.subplots(figsize=(10, 7))

    for file in files:

        df = pd.read_excel(file)
        df.columns = df.columns.map(str).str.strip()

        # Average displacement of columns 3 and 4 (metres)
        d3 = pd.to_numeric(df["3"], errors="coerce")
        d4 = pd.to_numeric(df["4"], errors="coerce")
        displacement = (d3 + d4) / 2

        # Zero displacement using mean of first 5 readings
        displacement = displacement - displacement.iloc[:5].mean()

        # Convert metres to millimetres and swap the sign
        displacement = -displacement * 1000

        # Calculate base shear (kN)
        base_shear = calculate_base_shear(df, heights)

        valid = displacement.notna() & base_shear.notna()
        x = displacement[valid]
        y = base_shear[valid]

        # Add each run to the specimen plot
        ax.plot(
            x, y,
            linewidth=0.7,
            label=file.stem.replace(f"Test{test_number}", "")
        )

    ax.set(
        xlabel="Displacement (mm)",
        ylabel="Base shear (kN)",
        title=f"{specimen} - All runs"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=3)

    fig.tight_layout()
    fig.savefig(
        OUTPUT / f"{specimen}_AllRuns.png",
        dpi=300
    )

    plt.show()
    plt.close(fig)

    print(f"Finished {specimen}: {len(files)} runs")

print("All plots saved in:", OUTPUT)
