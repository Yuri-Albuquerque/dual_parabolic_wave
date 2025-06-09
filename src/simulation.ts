import { Point2D, WaveParams, SimulationConfig } from './types';
import { Parabola } from './parabola';
import { WaveField } from './waveField';
import { OptimizedWaveField } from './optimizedWaveField';

export class DualParabolicWaveSimulation {
  private majorParabola!: Parabola;
  private minorParabola!: Parabola;
  private waveField!: OptimizedWaveField;
  private focusPoint!: Point2D;
  private config!: SimulationConfig;
  private waveParams!: WaveParams;
  private animationId: number | null = null;
  private canvas: HTMLCanvasElement | null = null;
  private ctx: CanvasRenderingContext2D | null = null;
  private lastRenderTime: number = 0;
  private renderInterval: number = 1000 / 30; // 30 FPS for rendering
  private simulationStepsPerFrame: number = 10; // Multiple physics steps per render

  constructor() {
    this.setupParabolas();
    this.setupWaveParameters();
    this.setupSimulationConfig();
    this.initializeWaveField();
  }

  private setupParabolas(): void {
    // Convert units: 20 inches = 508mm, focus distances in mm
    const majorDiameter = 20 * 25.4; // 508mm
    const majorFocus = 100; // mm
    const minorDiameter = 100; // mm
    const minorFocus = 50; // mm

    // Focus point is at origin (coincident focus)
    this.focusPoint = { x: 0, y: 0 };

    // Major parabola (umbrella, concave down)
    // Vertex is 100mm above focus point
    const majorVertex: Point2D = { x: 0, y: majorFocus };
    this.majorParabola = new Parabola(majorDiameter, majorFocus, majorVertex, false);

    // Minor parabola (bowl, concave up)
    // Vertex is 50mm below focus point
    const minorVertex: Point2D = { x: 0, y: -minorFocus };
    this.minorParabola = new Parabola(minorDiameter, minorFocus, minorVertex, true);
  }

  private setupWaveParameters(): void {
    const frequency = 1000; // Hz
    const speed = 343000; // mm/s (speed of sound in air)
    const wavelength = speed / frequency;

    this.waveParams = {
      frequency,
      wavelength,
      speed,
      amplitude: 1.0
    };
  }

  private setupSimulationConfig(): void {
    // Simulation domain should encompass both parabolas
    const xRange: [number, number] = [-300, 300];
    const yRange: [number, number] = [-100, 150];
    const gridSize = 300;

    // Calculate CFL-compliant time step
    const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
    const dy = (yRange[1] - yRange[0]) / (gridSize - 1);
    const minGridSpacing = Math.min(dx, dy);
    const speed = this.waveParams.speed;
    
    // CFL condition: dt <= CFL_factor * min(dx,dy) / (c * sqrt(2))
    // Use CFL_factor = 0.4 for stability margin
    const maxStableTimeStep = 0.4 * minGridSpacing / (speed * Math.sqrt(2));

    this.config = {
      gridSize,
      xRange,
      yRange,
      timeStep: maxStableTimeStep,
      dampingFactor: 0.001, // Minimal damping for better wave propagation
      reflectionCoeff: 0.95 // High reflection coefficient
    };

    console.log(`CFL-compliant time step: ${maxStableTimeStep.toExponential(2)} s`);
    console.log(`Grid spacing: dx=${dx.toFixed(3)}mm, dy=${dy.toFixed(3)}mm`);
  }

  private initializeWaveField(): void {
    this.waveField = new OptimizedWaveField(
      this.config,
      this.waveParams,
      this.focusPoint
    );
    
    // Set up update callback for asynchronous updates
    this.waveField.onUpdate((grid: Float32Array, time: number) => {
      // This will be called when the worker completes an update
      // No need to do anything here as the grid is automatically updated
    });
  }

  public setupVisualization(canvasId: string): void {
    this.canvas = document.getElementById(canvasId) as HTMLCanvasElement;
    if (!this.canvas) {
      throw new Error(`Canvas with id "${canvasId}" not found`);
    }

    this.ctx = this.canvas.getContext('2d');
    if (!this.ctx) {
      throw new Error('Could not get 2D context from canvas');
    }

    // Set canvas size
    this.canvas.width = 800;
    this.canvas.height = 600;
  }

  private drawParabolas(): void {
    if (!this.ctx || !this.canvas) return;

    const { xRange, yRange } = this.config;
    const canvasWidth = this.canvas.width;
    const canvasHeight = this.canvas.height;

    // Coordinate transformation functions
    const xToCanvas = (x: number) => ((x - xRange[0]) / (xRange[1] - xRange[0])) * canvasWidth;
    const yToCanvas = (y: number) => ((yRange[1] - y) / (yRange[1] - yRange[0])) * canvasHeight;

    this.ctx.strokeStyle = '#333';
    this.ctx.lineWidth = 2;

    // Draw major parabola
    this.ctx.beginPath();
    const majorParams = this.majorParabola.getParams();
    for (let x = -majorParams.diameter / 2; x <= majorParams.diameter / 2; x += 5) {
      const y = this.majorParabola.getY(x);
      const canvasX = xToCanvas(x);
      const canvasY = yToCanvas(y);
      
      if (x === -majorParams.diameter / 2) {
        this.ctx.moveTo(canvasX, canvasY);
      } else {
        this.ctx.lineTo(canvasX, canvasY);
      }
    }
    this.ctx.stroke();

    // Draw minor parabola
    this.ctx.beginPath();
    const minorParams = this.minorParabola.getParams();
    for (let x = -minorParams.diameter / 2; x <= minorParams.diameter / 2; x += 5) {
      const y = this.minorParabola.getY(x);
      const canvasX = xToCanvas(x);
      const canvasY = yToCanvas(y);
      
      if (x === -minorParams.diameter / 2) {
        this.ctx.moveTo(canvasX, canvasY);
      } else {
        this.ctx.lineTo(canvasX, canvasY);
      }
    }
    this.ctx.stroke();

    // Draw focus point
    this.ctx.fillStyle = '#ff0000';
    this.ctx.beginPath();
    this.ctx.arc(xToCanvas(this.focusPoint.x), yToCanvas(this.focusPoint.y), 5, 0, 2 * Math.PI);
    this.ctx.fill();
  }

  private drawWaveField(): void {
    if (!this.ctx || !this.canvas) return;

    const grid = this.waveField.getGrid();
    const boundaryMask = this.waveField.getBoundaryMask();
    const { gridSize } = this.config;
    const canvasWidth = this.canvas.width;
    const canvasHeight = this.canvas.height;

    const imageData = this.ctx.createImageData(canvasWidth, canvasHeight);
    const data = imageData.data;

    // Scale factors for mapping grid to canvas
    const xScale = canvasWidth / gridSize;
    const yScale = canvasHeight / gridSize;

    // Find min/max values for normalization
    let minVal = Infinity;
    let maxVal = -Infinity;
    for (let i = 0; i < gridSize; i++) {
      for (let j = 0; j < gridSize; j++) {
        const index = i * gridSize + j;
        if (boundaryMask[index]) {
          minVal = Math.min(minVal, grid[index]);
          maxVal = Math.max(maxVal, grid[index]);
        }
      }
    }

    const range = Math.max(Math.abs(minVal), Math.abs(maxVal));

    // Draw wave field
    for (let i = 0; i < gridSize; i++) {
      for (let j = 0; j < gridSize; j++) {
        const canvasX = Math.floor(j * xScale);
        const canvasY = Math.floor(i * yScale);

        if (canvasX >= canvasWidth || canvasY >= canvasHeight) continue;

        const pixelIndex = (canvasY * canvasWidth + canvasX) * 4;
        const gridIndex = i * gridSize + j;

        if (boundaryMask[gridIndex] && range > 0) {
          const normalizedValue = grid[gridIndex] / range;
          const intensity = Math.abs(normalizedValue);
          
          if (normalizedValue > 0) {
            // Positive values in red
            data[pixelIndex] = Math.floor(255 * intensity);     // R
            data[pixelIndex + 1] = 0;                           // G
            data[pixelIndex + 2] = 0;                           // B
          } else {
            // Negative values in blue
            data[pixelIndex] = 0;                               // R
            data[pixelIndex + 1] = 0;                           // G
            data[pixelIndex + 2] = Math.floor(255 * intensity); // B
          }
          data[pixelIndex + 3] = Math.floor(255 * intensity); // Alpha
        } else {
          // Outside boundary or zero amplitude
          data[pixelIndex] = 0;     // R
          data[pixelIndex + 1] = 0; // G
          data[pixelIndex + 2] = 0; // B
          data[pixelIndex + 3] = 0; // Alpha (transparent)
        }
      }
    }

    this.ctx.putImageData(imageData, 0, 0);
  }

  private render(): void {
    if (!this.ctx || !this.canvas) return;

    // Clear canvas
    this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

    // Draw wave field first (background)
    this.drawWaveField();

    // Draw parabolas on top
    this.drawParabolas();

    // Draw time and simulation info
    this.ctx.fillStyle = '#000';
    this.ctx.font = '16px Arial';
    this.ctx.fillText(`Time: ${this.waveField.getTime().toFixed(3)}s`, 10, 30);
    this.ctx.fillText(`Frequency: ${this.waveParams.frequency}Hz`, 10, 50);
  }

  private animate = (): void => {
    const currentTime = performance.now();
    
    // Control rendering frame rate independently of simulation
    if (currentTime - this.lastRenderTime >= this.renderInterval) {
      this.render();
      this.lastRenderTime = currentTime;
    }
    
    // Update simulation with multiple steps for efficiency
    if (this.waveField.isReady()) {
      this.waveField.update(this.config.timeStep, this.simulationStepsPerFrame);
    }
    
    this.animationId = requestAnimationFrame(this.animate);
  };

  public start(): void {
    if (this.animationId !== null) {
      this.stop();
    }
    this.animate();
  }

  public stop(): void {
    if (this.animationId !== null) {
      cancelAnimationFrame(this.animationId);
      this.animationId = null;
    }
  }

  public reset(): void {
    this.stop();
    this.waveField.reset();
    this.render();
  }

  public setFrequency(frequency: number): void {
    this.waveParams.frequency = frequency;
    this.waveParams.wavelength = this.waveParams.speed / frequency;
    this.reset();
  }

  public setAmplitude(amplitude: number): void {
    this.waveParams.amplitude = amplitude;
  }

  public destroy(): void {
    this.stop();
    this.waveField.destroy();
  }

  // Getters for accessing simulation parameters
  public getFocusPoint(): Point2D {
    return { ...this.focusPoint };
  }

  public getMajorParabola(): Parabola {
    return this.majorParabola;
  }

  public getMinorParabola(): Parabola {
    return this.minorParabola;
  }

  public getWaveParams(): WaveParams {
    return { ...this.waveParams };
  }
}
