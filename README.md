# simpulse
An single and period pulse injector for multiple purposes. Originally designed to simulated data to test FREDDA, by Harry Qiu (@hqiu-nju)

## Installation 
From source: 
```
git clone https://github.com/OwenJohnsons/simpulse.git
cd simpulse
pip install -e .
```

## Command-Line Tools 
### `simpulse` - Single Pulse Injector

Example: 
```
simpulse -o test --dm 500 --sig 3 --mode boxcar
```

Key options:
| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-o`, `--output` | str | `"test"` | Output file name prefix (no suffix). |
| `-m`, `--mode` | str | `"single"` | Pulse model: `single`, `scat`, `boxcar`. |
| `--snmode` | str | `"fluence"` | Injection mode: `fluence` or `snr`. |
| `-A`, `--amplitude` | float | `50` | Pulse amplitude (if not S/N scaling). |
| `-t`, `--tbin` | int | `10` | Time binning factor inside simulator. |
| `-f`, `--fbin` | int | `10` | Frequency binning factor inside simulator. |
| `-s`, `--samples` | int | `20000` | Number of samples per injected block. |
| `--nchan` | int | `336` | Number of frequency channels. |
| `--tsamp` | float | `1` | Sampling time in **ms**. |
| `--fch1` | float | `1100` | Center frequency of first channel (MHz). |
| `--bwchan` | float | `1` | Channel bandwidth (MHz). |
| `-N`, `--npulse` | int | `50` | Number of pulses inserted for each parameter combo. |
| `--dm_start` | float | `0` | Minimum DM value. |
| `--step` | float | `50` | DM step size. |
| `--dm` | float | `3000` | Maximum DM. |
| `--sig_start` | float | `0.5` | Minimum intrinsic sigma width (ms). |
| `--sig_step` | float | `0.5` | Width step size (ms). |
| `--sig` | float | `0.5` | Maximum intrinsic width (ms). |


### `simperiod` - Period Pulse Injector 

Example: 
```
simperiod -o pulsar --dm 50 --period 0.5 --pdot 1e-15 --width 5 --snr 100 --npulses 100
```

Key Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--dm` | float | *required* | Dispersion measure (pc cm⁻³). |
| `-p`, `--period` | float | *required* | Spin period **P** in seconds. |
| `--pdot` | float | `0.0` | Period derivative **Pdot** (s/s). |
| `-w`, `--width` | float | *required* | Intrinsic Gaussian pulse width (ms). |
| `--snr` | float | *required* | Target integrated S/N (L2\_clean scaling). |
| `--npulses` | int | `100` | Number of pulses to simulate. |
| `--fch1` | float | `1100` | Top-of-band channel center frequency (MHz). |
| `--bwchan` | float | `1.0` | Channel bandwidth (MHz). |
| `--nchan` | int | `336` | Number of frequency channels. |
| `--tsamp` | float | `1.0` | Sampling time in **ms**. |
| `--tbin` | int | `10` | Internal time binning. |
| `--fbin` | int | `10` | Internal frequency binning. |
| `--noise-std` | float | `18.0` | Noise standard deviation. |
| `--noise-base` | float | `127.0` | Noise baseline before uint8 casting. |
| `-o`, `--outfile` | str | `"simperiodic"` | Output `.fil` filename (without extension). |


## Authors 
Harry Qiu (SKAO), original author 

Owen A. Johnson (TCD), housekeeping 
