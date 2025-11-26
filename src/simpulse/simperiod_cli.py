'''
Code Purpose: Generate synthetic periodic signals for testing search pipelines and calibrtation workflows.
Author: Owen A. Johnson 
Date: November 2025
'''
# external imports
import argparse
import numpy as np
from rich.console import Console
from rich.progress import Progress, BarColumn, TextColumn, TimeRemainingColumn

# internal imports
from simpulse.sim.model import Spectra
from simpulse.sim.burst import tidm
from simpulse.sim.measurement import L2_clean
from simpulse.io.fbio import makefilterbank

console = Console()

def main():
    parser = argparse.ArgumentParser(
        description="Simulate a periodic pulsar with given DM, period, Pdot, width, and target S/N."
    )

    # Pulsar parameters
    parser.add_argument("-dm", type=float, required=True,
                        help="Dispersion measure (pc cm^-3)")
    parser.add_argument("-p", "--period", type=float, required=True,
                        help="Spin period P (seconds)")
    parser.add_argument("-pdot", type=float, default=0.0,
                        help="Period derivative Pdot (s/s)")
    parser.add_argument("-w", "--width", type=float, required=True,
                        help="Intrinsic pulse width sigma (ms)")
    parser.add_argument("-snr", type=float, required=True,
                        help="Target integrated S/N (L2_clean definition)")
    parser.add_argument("-npulses", type=int, default=100,
                        help="Number of pulses to simulate")

    # Data parameters
    parser.add_argument("-fch1", type=float, default=190,
                        help="Top-of-band channel centre frequency (MHz)")
    parser.add_argument("-bwchan", type=float, default=-0.02435,
                        help="Channel bandwidth (MHz)")
    parser.add_argument("-nchan", type=int, default=3296,
                        help="Number of frequency channels")
    parser.add_argument("-tsamp", type=float, default=0.655,
                        help="Time resolution (ms)")
    parser.add_argument("-tbin", type=int, default=10,
                        help="Time grid resolution (used internally)")
    parser.add_argument("-fbin", type=int, default=10,
                        help="Frequency grid resolution (used internally)")

    # Noise + output
    parser.add_argument("--noise-std", type=float, default=18.0,
                        help="Standard deviation of background noise (pre-scaling)")
    parser.add_argument("--noise-base", type=float, default=127.0,
                        help="Base level added to data before uint8 cast")
    parser.add_argument("-o", "--output", type=str, default="simperiodic",
                        help="Output filterbank basename ('.fil' will be added)")

    args = parser.parse_args()

    simulate_periodic(args)
    
def simulate_periodic(args):
    console.print("[bold magenta]SIMPERIOD pulsar simulator[/]")

    dm = args.dm
    P0_s = args.period
    pdot_s = args.pdot
    width_ms = args.width
    target_snr = args.snr
    npulses = args.npulses

    tsamp_ms = args.tsamp
    fch1 = args.fch1
    bwchan = args.bwchan
    nchan = args.nchan
    tbin = args.tbin
    fbin = args.fbin

    noise_std = args.noise_std
    noise_base = args.noise_base
    output = args.output

    # --- 1. Compute pulse emission times including Pdot ---
    # t_n = n*P0 + 0.5*n*(n-1)*Pdot  (seconds)
    console.print(f"[bold]DM[/]: {dm} pc cm^-3")
    console.print(f"[bold]P[/]: {P0_s} s, [bold]Pdot[/]: {pdot_s} s/s")
    console.print(f"[bold]Width[/]: {width_ms} ms | [bold]Pulses[/]: {npulses}")
    console.print(f"[bold]Target S/N[/]: {target_snr}\n")

    n_arr = np.arange(npulses, dtype=float)
    t_n_s = n_arr * P0_s + 0.5 * n_arr * (n_arr - 1.0) * pdot_s
    t_last_s = t_n_s[-1]
    total_s = t_last_s + 2.0 * P0_s  # small padding
    total_ms = total_s * 1000.0

    nsamp = int(total_ms / tsamp_ms) + 1
    console.print(f"Total duration: [bold]{total_s:.2f} s[/] "
                  f"({nsamp} samples at {tsamp_ms} ms)\n")

    # --- 2. Set up Spectra just to get header + freq grid ---
    spec = Spectra(fch1=fch1, nchan=nchan, bwchan=bwchan,
                   tsamp=tsamp_ms, tbin=tbin, fbin=fbin)

    vif = spec.vif  # frequency grid (MHz)
    nchan = spec.nchan

    # --- 3. Build burst-only dynamic spectrum (no noise yet) ---
    t = np.arange(nsamp) * tsamp_ms  # time axis in ms
    burst_dyn = np.zeros((nsamp, nchan), dtype=float)

    console.print("[bold blue]Injecting pulses...[/]")

    with Progress(
        TextColumn("[cyan]{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total} pulses"),
        TimeRemainingColumn(),
        console=console,
    ) as progress:

        pulse_task = progress.add_task("Pulses", total=npulses)

        for n_idx, t_emit_s in enumerate(t_n_s):
            t_emit_ms = t_emit_s * 1000.0

            # For each channel, compute arrival time including DM delay
            for ci, freq in enumerate(vif):
                t_arrival_ms = t_emit_ms + tidm(dm, freq, fch1)
                burst_dyn[:, ci] += np.exp(
                    -0.5 * ((t - t_arrival_ms) / width_ms) ** 2
                )

            progress.update(pulse_task, advance=1)

    # --- 4. Scale burst to target S/N using L2_clean ---
    console.print("\n[bold blue]Measuring clean S/N for scaling...[/]")
    snr0 = L2_clean(burst_dyn)
    if snr0 == 0:
        console.print("[bold red]Error:[/] clean S/N is zero; check parameters.")
        return

    amp_factor = target_snr / snr0
    console.print(f"Base clean S/N: {snr0:.2f} → scaling by factor {amp_factor:.3f}")
    burst_dyn *= amp_factor

    # --- 5. Add noise + base level ---
    console.print("[bold blue]Adding noise and base level...[/]")
    noise = np.random.randn(nsamp, nchan) * noise_std + noise_base
    dyn = noise + burst_dyn

    # --- 6. Write to filterbank ---
    console.print(f"[bold blue]Writing filterbank:[/] {output}.fil")

    header = spec.header.copy()
    header["nsamples"] = nsamp

    fbank = makefilterbank(output + ".fil", header=header)
    fbank.writeblock(dyn.astype(np.uint8))
    fbank.closefile()

    console.print("[bold green]Done![/] Filterbank written.\n")