import { DualParabolicWaveSimulation } from './simulation';

// For Node.js environment, we'll create a simple console-based visualization
export class ConsoleRenderer {
  private simulation: DualParabolicWaveSimulation;
  private intervalId: NodeJS.Timeout | null = null;

  constructor() {
    this.simulation = new DualParabolicWaveSimulation();
  }

  public renderToConsole(): void {
    const grid = this.simulation['waveField'].getGrid();
    const boundaryMask = this.simulation['waveField'].getBoundaryMask();
    const gridSize = this.simulation['config'].gridSize;
    
    // Clear console
    console.clear();
    
    console.log('='.repeat(60));
    console.log('Dual Parabolic Wave Propagation Simulation');
    console.log(`Time: ${this.simulation['waveField'].getTime().toFixed(3)}s`);
    console.log(`Frequency: ${this.simulation.getWaveParams().frequency}Hz`);
    console.log('='.repeat(60));

    // Sample every few grid points for console display
    const step = Math.floor(gridSize / 40);
    let output = '';

    for (let i = 0; i < gridSize; i += step) {
      let line = '';
      for (let j = 0; j < gridSize; j += step) {
        const index = i * gridSize + j;
        if (boundaryMask[index]) {
          const value = grid[index];
          if (Math.abs(value) < 0.1) {
            line += ' ';
          } else if (value > 0) {
            line += '+';
          } else {
            line += '-';
          }
        } else {
          line += '█'; // Boundary
        }
      }
      output += line + '\n';
    }
    
    console.log(output);
    
    // Print some statistics
    let maxAmplitude = 0;
    let activePoints = 0;
    
    for (let i = 0; i < gridSize; i++) {
      for (let j = 0; j < gridSize; j++) {
        const index = i * gridSize + j;
        if (boundaryMask[index]) {
          activePoints++;
          maxAmplitude = Math.max(maxAmplitude, Math.abs(grid[index]));
        }
      }
    }
    
    console.log(`Active grid points: ${activePoints}`);
    console.log(`Max amplitude: ${maxAmplitude.toFixed(4)}`);
    console.log('Legend: █=boundary, +=positive wave, -=negative wave, space=neutral');
  }

  public startConsoleAnimation(intervalMs: number = 100): void {
    this.intervalId = setInterval(() => {
      this.simulation['waveField'].update(0.001);
      this.renderToConsole();
    }, intervalMs);
  }

  public stopConsoleAnimation(): void {
    if (this.intervalId) {
      clearInterval(this.intervalId);
      this.intervalId = null;
    }
  }

  public getSimulation(): DualParabolicWaveSimulation {
    return this.simulation;
  }
}

// Export everything for use in browser or Node.js
export { DualParabolicWaveSimulation } from './simulation';
export { Parabola } from './parabola';
export { WaveField } from './waveField';
export * from './types';

// Main execution for Node.js
if (require.main === module) {
  console.log('Starting Dual Parabolic Wave Simulation...\n');
  console.log('Simulation Parameters:');
  console.log('- Major parabola: 20 inches (508mm) diameter, 100mm focus');
  console.log('- Minor parabola: 100mm diameter, 50mm focus');
  console.log('- Wave frequency: 1000Hz');
  console.log('- Sound speed: 343 m/s');
  console.log('\nPress Ctrl+C to stop the simulation\n');

  const renderer = new ConsoleRenderer();
  renderer.startConsoleAnimation(200);

  // Graceful shutdown
  process.on('SIGINT', () => {
    console.log('\nStopping simulation...');
    renderer.stopConsoleAnimation();
    process.exit(0);
  });
}
