# Dual Parabolic Wave Propagation Simulation

This project implements a wave propagation simulation inside dual parabolic constraints using TypeScript.

## System Description

The simulation models an acoustic system with two parabolic reflectors:

- **Major Parabola (Umbrella)**: 20 inches (508mm) diameter, concave down, focus 100mm from vertex
- **Minor Parabola (Bowl)**: 100mm diameter, concave up, focus 50mm from vertex
- **Coincident Focus**: Both parabolas share the same focus point where the sound source is placed
- **Wave Source**: Point source generating sinusoidal waves at the focus

## Installation and Setup

### Prerequisites

1. **Node.js and npm**: Install from [nodejs.org](https://nodejs.org/)
   ```bash
   # On Ubuntu/Debian:
   sudo apt update
   sudo apt install nodejs npm
   
   # On macOS with Homebrew:
   brew install node
   
   # On Windows: Download installer from nodejs.org
   ```

2. **TypeScript**: Will be installed as a dependency

### Installation Steps

1. **Clone or navigate to the project directory**:
   ```bash
   cd /path/to/dual_parabolic_wave
   ```

2. **Install dependencies**:
   ```bash
   npm install
   ```

3. **Build the TypeScript code**:
   ```bash
   npm run build
   ```

4. **Run the simulation**:
   
   **Option A: Console version (Node.js)**
   ```bash
   npm start
   # or
   npm run dev
   ```
   
   **Option B: Browser version**
   Open `index.html` in a web browser after building the project.

## Usage

### Console Version
The console version provides a text-based visualization of the wave propagation:
- `+` represents positive wave amplitude
- `-` represents negative wave amplitude
- `█` represents parabolic boundaries
- Space represents neutral/zero amplitude

### Browser Version
The browser version provides a real-time graphical visualization with:
- Red areas: Positive wave amplitude
- Blue areas: Negative wave amplitude
- Black lines: Parabolic boundaries
- Red dot: Sound source at focus point

### Controls
- **Start/Stop**: Control animation
- **Reset**: Reset simulation to initial state
- **Frequency**: Adjust wave frequency (100-3000 Hz)
- **Amplitude**: Adjust wave amplitude (0.1-2.0)

## Project Structure

```
src/
├── types.ts          # TypeScript interfaces and type definitions
├── parabola.ts       # Parabola geometry and reflection calculations
├── waveField.ts      # Wave equation solver and field management
├── simulation.ts     # Main simulation class with visualization
└── index.ts          # Entry point and console renderer

index.html            # Browser-based visualization interface
package.json          # Node.js dependencies and scripts
tsconfig.json         # TypeScript configuration
```

## Physics Implementation

### Wave Equation
The simulation solves the 2D wave equation:
```
∂²u/∂t² = c²∇²u
```
Where:
- `u` is the wave field
- `c` is the speed of sound (343 m/s)
- `∇²` is the Laplacian operator

### Parabolic Geometry
- **Major parabola**: `y = -x²/(4f) + v_y` (concave down)
- **Minor parabola**: `y = x²/(4f) + v_y` (concave up)
- Where `f` is the focal length and `v_y` is the vertex y-coordinate

### Boundary Conditions
- **Rigid boundaries**: Zero displacement at parabolic surfaces
- **Reflection**: Acoustic waves reflect from parabolic surfaces
- **Absorption**: Small damping factor to prevent numerical instabilities

### Numerical Methods
- **Finite Difference**: 2nd order accurate in space and time
- **Grid Resolution**: Configurable grid size (default 400x400)
- **Time Stepping**: Explicit time integration with stability constraints

## Customization

### Modify Parabola Parameters
Edit the `setupParabolas()` method in `src/simulation.ts`:
```typescript
const majorDiameter = 20 * 25.4; // Change diameter
const majorFocus = 100;           // Change focus distance
```

### Adjust Wave Parameters
Edit the `setupWaveParameters()` method:
```typescript
const frequency = 1000;  // Change frequency (Hz)
const speed = 343000;    // Change sound speed (mm/s)
```

### Simulation Settings
Edit the `setupSimulationConfig()` method:
```typescript
gridSize: 400,           // Grid resolution
timeStep: 0.001,         // Time step (seconds)
dampingFactor: 0.001,    // Damping coefficient
reflectionCoeff: 0.95    // Reflection coefficient
```

## Troubleshooting

### Common Issues

1. **"npm: command not found"**
   - Install Node.js and npm first
   - Verify installation: `node --version && npm --version`

2. **TypeScript compilation errors**
   - Ensure TypeScript is installed: `npm install -g typescript`
   - Check tsconfig.json settings

3. **Canvas not displaying in browser**
   - Ensure the project is built: `npm run build`
   - Check browser console for JavaScript errors
   - Try serving the file through a local server

4. **Performance issues**
   - Reduce grid size in simulation config
   - Increase time step (less accurate but faster)
   - Lower animation frame rate

### Performance Optimization

For better performance:
- Reduce `gridSize` to 200 or 300
- Increase `timeStep` to 0.002 or higher
- Add more damping to reduce computational load
- Use Web Workers for background computation (advanced)

## Scientific Applications

This simulation can be used to study:
- **Acoustic focusing**: How parabolic reflectors concentrate sound waves
- **Wave interference**: Constructive and destructive interference patterns
- **Resonance**: Standing wave patterns in parabolic cavities
- **Acoustic design**: Optimization of parabolic reflector systems

## License

[Your chosen license - MIT recommended for open source]
This code implements a dual parabolic wave propagation in TypeScript
