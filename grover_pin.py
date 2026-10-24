"""Problem 17: find the hidden 4-bit PIN 1011 with Grover's algorithm (Qiskit + Aer)."""
import math

import matplotlib
matplotlib.use("Agg")  # write PNGs without needing a display
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, transpile
from qiskit.visualization import plot_histogram
from qiskit_aer import AerSimulator

PIN = "1011"
N_QUBITS = len(PIN)
SHOTS = 1024
SIM = AerSimulator()


def oracle():
    """Phase-flip only |PIN>: |x> -> -|x> if x == PIN, else unchanged."""
    qc = QuantumCircuit(N_QUBITS, name="Oracle")
    # Qiskit is little-endian: qubit 0 is the RIGHTMOST bit of a bitstring.
    # Reversing PIN lines up character i with qubit i.
    zeros = [i for i, bit in enumerate(reversed(PIN)) if bit == "0"]
    if zeros:
        qc.x(zeros)  # map |PIN> to |1111>
    # Multi-controlled Z = H . MCX . H on the last qubit; flips the phase of |1111>.
    qc.h(N_QUBITS - 1)
    qc.mcx(list(range(N_QUBITS - 1)), N_QUBITS - 1)
    qc.h(N_QUBITS - 1)
    if zeros:
        qc.x(zeros)  # undo the mapping
    return qc


def diffuser():
    """Inversion about the mean: 2|s><s| - I, where |s> is the uniform superposition."""
    qc = QuantumCircuit(N_QUBITS, name="Diffuser")
    qc.h(range(N_QUBITS))
    qc.x(range(N_QUBITS))  # |0000> -> |1111>
    qc.h(N_QUBITS - 1)
    qc.mcx(list(range(N_QUBITS - 1)), N_QUBITS - 1)  # phase-flip |1111> (i.e. |0000>)
    qc.h(N_QUBITS - 1)
    qc.x(range(N_QUBITS))
    qc.h(range(N_QUBITS))
    return qc


def grover(k):
    """Uniform superposition, then k rounds of (oracle, diffuser), then measure."""
    qc = QuantumCircuit(N_QUBITS, N_QUBITS)
    qc.h(range(N_QUBITS))
    for _ in range(k):
        qc.compose(oracle(), inplace=True)
        qc.compose(diffuser(), inplace=True)
    qc.measure(range(N_QUBITS), range(N_QUBITS))  # qubit i -> classical bit i
    return qc


def run(k):
    """Return measurement counts for k Grover iterations."""
    job = SIM.run(transpile(grover(k), SIM), shots=SHOTS, seed_simulator=42)
    return job.result().get_counts()


def theory(k):
    """Ideal success probability sin^2((2k+1) * theta), with sin(theta) = 1/sqrt(N) = 1/4."""
    theta = math.asin(1 / math.sqrt(2 ** N_QUBITS))
    return math.sin((2 * k + 1) * theta) ** 2


def main():
    # --- Part 1: 3 iterations (optimal: floor(pi/4 * sqrt(16)) = 3) ---
    counts = run(3)
    top = max(counts, key=counts.get)
    print(f"k=3 counts: {dict(sorted(counts.items()))}")
    print(f"Most frequent result: {top}  ({counts[top]}/{SHOTS})")
    assert top == PIN, "Grover did not find the PIN - check oracle bit order"

    fig = plot_histogram(counts, title=f"Grover search for PIN {PIN} (3 iterations, {SHOTS} shots)",
                         figsize=(10, 5))
    fig.savefig("grover_histogram.png", dpi=150, bbox_inches="tight")

    # --- Part 2: sweep k = 1..6 and compare with theory ---
    ks = list(range(1, 7))
    sim_p = [run(k).get(PIN, 0) / SHOTS for k in ks]
    th_p = [theory(k) for k in ks]

    print("\n k | simulated | theory | diff")
    for k, s, t in zip(ks, sim_p, th_p):
        print(f" {k} |   {s:.4f}  | {t:.4f} | {s - t:+.4f}")
        # 1024 shots -> std error <= 0.016; 0.05 is ~3 sigma
        assert abs(s - t) < 0.05, f"k={k}: simulation disagrees with theory"
    assert max(ks, key=lambda k: sim_p[k - 1]) == 3, "Peak should be at k=3"

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(ks, sim_p, "o-", label=f"Simulated (Aer, {SHOTS} shots)")
    ax.plot(ks, th_p, "s--", label=r"Theory $\sin^2((2k+1)\theta)$, $\theta=\arcsin(1/4)$")
    ax.axvline(3, color="gray", linestyle=":", label="Optimal k = 3")
    ax.set(xlabel="Grover iterations k", ylabel=f"P(measure {PIN})", ylim=(0, 1.05),
           title=f"Success probability vs iterations (PIN {PIN})")
    ax.set_xticks(ks)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig("grover_success_vs_iterations.png", dpi=150, bbox_inches="tight")

    print("\nSaved grover_histogram.png and grover_success_vs_iterations.png")


if __name__ == "__main__":
    main()
