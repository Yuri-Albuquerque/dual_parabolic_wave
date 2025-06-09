import { Point2D, WaveParams, SimulationConfig } from './types';
import { Parabola } from './parabola';

export class WaveField {
  private grid!: number[][];
  private previousGrid!: number[][];
  private velocityGrid!: number[][];
  private boundaryMask!: boolean[][];
  private config: SimulationConfig;
  private waveParams: WaveParams;
  private majorParabola: Parabola;
  private minorParabola: Parabola;
  private focusPoint: Point2D;
  private time: number;

  constructor(
    config: SimulationConfig,
    waveParams: WaveParams,
    majorParabola: Parabola,
    minorParabola: Parabola,
    focusPoint: Point2D
  ) {
    this.config = config;
    this.waveParams = waveParams;
    this.majorParabola = majorParabola;
    this.minorParabola = minorParabola;
    this.focusPoint = focusPoint;
    this.time = 0;

    this.validateCFLCondition();
    this.initializeGrids();
    this.createBoundaryMask();
  }

  private validateCFLCondition(): void {
    const { gridSize, xRange, yRange, timeStep } = this.config;
    const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
    const dy = (yRange[1] - yRange[0]) / (gridSize - 1);
    const minGridSpacing = Math.min(dx, dy);
    
    // CFL condition for 2D wave equation: Δt ≤ CFL_factor × min(Δx, Δy) / (c × √2)
    // We use CFL_factor = 0.5 for stability
    const maxStableTimeStep = 0.5 * minGridSpacing / (this.waveParams.speed * Math.sqrt(2));
    
    if (timeStep > maxStableTimeStep) {
      console.warn(`Time step ${timeStep} exceeds CFL stability limit ${maxStableTimeStep}`);
      console.warn(`Grid spacing: dx=${dx}, dy=${dy}, wave speed=${this.waveParams.speed}`);
      console.warn(`Recommended max time step: ${maxStableTimeStep}`);
    }
  }

  private initializeGrids(): void {
    const { gridSize } = this.config;
    
    this.grid = Array(gridSize).fill(null).map(() => Array(gridSize).fill(0));
    this.previousGrid = Array(gridSize).fill(null).map(() => Array(gridSize).fill(0));
    this.velocityGrid = Array(gridSize).fill(null).map(() => Array(gridSize).fill(0));
    this.boundaryMask = Array(gridSize).fill(null).map(() => Array(gridSize).fill(false));
  }

  private createBoundaryMask(): void {
    const { gridSize, xRange, yRange } = this.config;
    const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
    const dy = (yRange[1] - yRange[0]) / (gridSize - 1);

    for (let i = 0; i < gridSize; i++) {
      for (let j = 0; j < gridSize; j++) {
        const x = xRange[0] + j * dx;
        const y = yRange[1] - i * dy; // Flip y-axis for proper visualization

        const point: Point2D = { x, y };
        
        // Point is inside the cavity if it's inside the major parabola AND outside the minor parabola
        const insideMajor = this.majorParabola.isPointInside(point);
        const outsideMinor = !this.minorParabola.isPointInside(point);
        
        this.boundaryMask[i][j] = insideMajor && outsideMinor;
      }
    }
  }

  private getGridCoordinates(gridI: number, gridJ: number): Point2D {
    const { gridSize, xRange, yRange } = this.config;
    const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
    const dy = (yRange[1] - yRange[0]) / (gridSize - 1);

    return {
      x: xRange[0] + gridJ * dx,
      y: yRange[1] - gridI * dy
    };
  }

  private calculateSourceAmplitude(point: Point2D, time: number): number {
    const distance = Math.sqrt(
      Math.pow(point.x - this.focusPoint.x, 2) + 
      Math.pow(point.y - this.focusPoint.y, 2)
    );

    if (distance < 1e-6) {
      return this.waveParams.amplitude * Math.sin(2 * Math.PI * this.waveParams.frequency * time);
    }

    const k = 2 * Math.PI / this.waveParams.wavelength;
    const phase = k * distance - 2 * Math.PI * this.waveParams.frequency * time;
    
    // Spherical wave with 1/r amplitude decay
    return (this.waveParams.amplitude / (distance + 1)) * Math.sin(phase);
  }

  private applyWaveEquation(dt: number): void {
    const { gridSize, dampingFactor, xRange, yRange } = this.config;
    const c = this.waveParams.speed; // wave speed
    const c2 = c * c;
    
    // Calculate grid spacing
    const hx = (xRange[1] - xRange[0]) / (gridSize - 1);
    const hz = (yRange[1] - yRange[0]) / (gridSize - 1);
    const ht = dt;

    // Pre-calculate coefficients (following C++ implementation)
    const q0 = c2 * ht;
    const q1 = c2 * ht * ht;
    const q2 = ((c * ht) / hx) * ((c * ht) / hx);
    const q3 = ((c * ht) / hz) * ((c * ht) / hz);

    // Finite difference coefficients for Laplacian
    const Cxx = [1.0, -2.0, 1.0];
    const Czz = [1.0, -2.0, 1.0];

    // Create new grid to avoid in-place updates
    const newGrid = Array(gridSize).fill(null).map(() => Array(gridSize).fill(0));

    // Wave equation implementation following C++ solver
    // u(t+1) = (-1.0 * (u(t-1) - 2.0 * u(t)) + q0 * damp * u(t-1) + q1 * source + q2 * dxx + q3 * dzz) / (1.0 + damp * q0)

    for (let i = 1; i < gridSize - 1; i++) {
      for (let j = 1; j < gridSize - 1; j++) {
        if (!this.boundaryMask[i][j]) {
          newGrid[i][j] = 0;
          continue;
        }

        // Calculate second derivatives using finite difference stencil
        const dxx = Cxx[0] * this.grid[i][j-1] + Cxx[1] * this.grid[i][j] + Cxx[2] * this.grid[i][j+1];
        const dzz = Czz[0] * this.grid[i-1][j] + Czz[1] * this.grid[i][j] + Czz[2] * this.grid[i+1][j];

        // Source term (will be added separately)
        const sourceValue = 0; // This is handled in addSourceExcitation

        // Apply wave equation with damping (following C++ implementation)
        const numerator = -1.0 * (this.previousGrid[i][j] - 2.0 * this.grid[i][j]) + 
                         q0 * dampingFactor * this.previousGrid[i][j] + 
                         q1 * sourceValue + 
                         q2 * dxx + 
                         q3 * dzz;
        
        const denominator = 1.0 + dampingFactor * q0;
        
        newGrid[i][j] = numerator / denominator;
      }
    }

    // Handle boundary conditions for edges (following C++ special case for i=0)
    for (let j = 1; j < gridSize - 1; j++) {
      if (this.boundaryMask[0][j]) {
        const dxx = Cxx[0] * this.grid[0][j-1] + Cxx[1] * this.grid[0][j] + Cxx[2] * this.grid[0][j+1];
        // Special boundary condition for top edge (i=0)
        const dzz_boundary = 2.0 * (this.grid[1][j] - this.grid[0][j]);
        
        const numerator = -1.0 * (this.previousGrid[0][j] - 2.0 * this.grid[0][j]) + 
                         q0 * dampingFactor * this.previousGrid[0][j] + 
                         q2 * dxx + 
                         q3 * dzz_boundary;
        
        const denominator = 1.0 + dampingFactor * q0;
        newGrid[0][j] = numerator / denominator;
      }
    }

    // Update grids
    this.previousGrid = this.grid.map(row => [...row]);
    this.grid = newGrid;

    // Apply boundary conditions
    this.applyBoundaryConditions();
  }

  private applyBoundaryConditions(): void {
    const { gridSize, reflectionCoeff } = this.config;

    for (let i = 0; i < gridSize; i++) {
      for (let j = 0; j < gridSize; j++) {
        if (!this.boundaryMask[i][j]) {
          this.grid[i][j] = 0;
          this.previousGrid[i][j] = 0;
          this.velocityGrid[i][j] = 0;
        } else {
          // Check for boundary neighbors and apply reflection
          const neighbors = [
            [i-1, j], [i+1, j], [i, j-1], [i, j+1]
          ];

          let isBoundaryPoint = false;
          for (const [ni, nj] of neighbors) {
            if (ni >= 0 && ni < gridSize && nj >= 0 && nj < gridSize) {
              if (!this.boundaryMask[ni][nj]) {
                isBoundaryPoint = true;
                break;
              }
            }
          }

          if (isBoundaryPoint) {
            this.grid[i][j] *= reflectionCoeff;
          }
        }
      }
    }
  }

  private addSourceExcitation(time: number): void {
    const { gridSize, xRange, yRange } = this.config;
    const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
    const dy = (yRange[1] - yRange[0]) / (gridSize - 1);

    // Find grid point closest to focus
    const focusI = Math.round((yRange[1] - this.focusPoint.y) / dy);
    const focusJ = Math.round((this.focusPoint.x - xRange[0]) / dx);

    if (focusI >= 0 && focusI < gridSize && focusJ >= 0 && focusJ < gridSize) {
      if (this.boundaryMask[focusI][focusJ]) {
        // Single pulse at initial time - Gaussian envelope with short duration
        const pulseWidth = 1.0 / this.waveParams.frequency; // One period width
        const pulseDuration = 2.0 * pulseWidth; // Duration of the pulse
        
        let sourceValue = 0.0;
        if (time <= pulseDuration) {
          // Gaussian envelope for smooth pulse
          const gaussianWidth = pulseWidth / 3.0;
          const envelope = Math.exp(-Math.pow(time - pulseWidth, 2) / (2 * gaussianWidth * gaussianWidth));
          
          // Single frequency pulse
          const sourceAmplitude = this.waveParams.amplitude * 10.0; // Stronger source
          sourceValue = sourceAmplitude * envelope * Math.sin(2 * Math.PI * this.waveParams.frequency * time);
        }
        
        // Apply source directly to the grid (this will be used in the next wave equation step)
        // The C++ implementation adds q1 * source to the numerator
        this.grid[focusI][focusJ] += sourceValue * 0.1; // Scale down to prevent instability
        
        // Add to neighboring points for better source distribution
        const neighbors = [
          [focusI-1, focusJ], [focusI+1, focusJ], 
          [focusI, focusJ-1], [focusI, focusJ+1]
        ];
        
        for (const [ni, nj] of neighbors) {
          if (ni >= 0 && ni < gridSize && nj >= 0 && nj < gridSize) {
            if (this.boundaryMask[ni][nj]) {
              this.grid[ni][nj] += sourceValue * 0.05; // Half strength for neighbors
            }
          }
        }
      }
    }
  }

  public update(dt: number): void {
    this.time += dt;
    
    // Add source excitation
    this.addSourceExcitation(this.time);
    
    // Apply wave equation
    this.applyWaveEquation(dt);
  }

  public getGrid(): number[][] {
    return this.grid.map(row => [...row]);
  }

  public getBoundaryMask(): boolean[][] {
    return this.boundaryMask.map(row => [...row]);
  }

  public getTime(): number {
    return this.time;
  }

  public reset(): void {
    this.time = 0;
    this.initializeGrids();
  }
}
