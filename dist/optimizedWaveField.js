"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.OptimizedWaveField = void 0;
class OptimizedWaveField {
    constructor(config, waveParams, focusPoint) {
        this.worker = null;
        this.time = 0;
        this.isInitialized = false;
        this.updateCallbacks = [];
        this.config = config;
        this.waveParams = waveParams;
        this.focusPoint = focusPoint;
        const gridSize = config.gridSize;
        this.grid = new Float32Array(gridSize * gridSize);
        this.boundaryMask = new Uint8Array(gridSize * gridSize);
        this.initializeWorker();
    }
    async initializeWorker() {
        try {
            // Create worker from the TypeScript file (will be compiled)
            const workerBlob = await this.createWorkerBlob();
            this.worker = new Worker(URL.createObjectURL(workerBlob));
            this.worker.onmessage = (event) => {
                const { type, data } = event.data;
                switch (type) {
                    case 'initialized':
                        this.grid = new Float32Array(data.grid);
                        this.boundaryMask = new Uint8Array(data.boundaryMask);
                        this.time = data.time;
                        this.isInitialized = true;
                        break;
                    case 'updated':
                        this.grid = new Float32Array(data.grid);
                        this.time = data.time;
                        // Notify callbacks
                        this.updateCallbacks.forEach(callback => callback(this.grid, this.time));
                        break;
                    case 'reset':
                        this.grid = new Float32Array(data.grid);
                        this.time = data.time;
                        break;
                    case 'error':
                        console.error('Worker error:', data.message);
                        break;
                }
            };
            this.worker.onerror = (error) => {
                console.error('Worker error:', error);
            };
            // Initialize the worker
            this.worker.postMessage({
                type: 'init',
                data: {
                    config: this.config,
                    waveParams: this.waveParams,
                    focusPoint: this.focusPoint
                }
            });
        }
        catch (error) {
            console.error('Failed to initialize worker, falling back to main thread');
            this.worker = null;
        }
    }
    async createWorkerBlob() {
        // Read the worker TypeScript file and compile it to JavaScript
        const workerCode = `
      // Inline worker code (simplified version)
      class WorkerWaveField {
        constructor(config, waveParams, focusPoint) {
          this.config = config;
          this.waveParams = waveParams;
          this.focusPoint = focusPoint;
          this.time = 0;
          
          const gridSize = config.gridSize;
          this.grid = new Float32Array(gridSize * gridSize);
          this.previousGrid = new Float32Array(gridSize * gridSize);
          this.boundaryMask = new Uint8Array(gridSize * gridSize);
          this.sourceGrid = new Float32Array(gridSize * gridSize);
          
          this.initializeBoundaryMask();
          this.calculateFocusPosition();
        }

        initializeBoundaryMask() {
          const { gridSize, xRange, yRange } = this.config;
          const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
          const dy = (yRange[1] - yRange[0]) / (gridSize - 1);

          const majorDiameter = 20 * 25.4;
          const majorFocus = 100;
          const minorDiameter = 100;
          const minorFocus = 50;

          for (let i = 0; i < gridSize; i++) {
            for (let j = 0; j < gridSize; j++) {
              const x = xRange[0] + j * dx;
              const y = yRange[1] - i * dy;

              const majorY = -(x * x) / (4 * majorFocus) + majorFocus;
              const insideMajor = y <= majorY && Math.abs(x) <= majorDiameter / 2;

              const minorY = (x * x) / (4 * minorFocus) - minorFocus;
              const outsideMinor = y >= minorY || Math.abs(x) > minorDiameter / 2;

              const index = i * gridSize + j;
              this.boundaryMask[index] = (insideMajor && outsideMinor) ? 1 : 0;
            }
          }
        }

        calculateFocusPosition() {
          const { gridSize, xRange, yRange } = this.config;
          const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
          const dy = (yRange[1] - yRange[0]) / (gridSize - 1);

          this.focusI = Math.round((yRange[1] - this.focusPoint.y) / dy);
          this.focusJ = Math.round((this.focusPoint.x - xRange[0]) / dx);
        }

        applyWaveEquation(dt) {
          const { gridSize, dampingFactor, xRange, yRange } = this.config;
          const c = this.waveParams.speed;
          const c2 = c * c;
          const dx = (xRange[1] - xRange[0]) / (gridSize - 1);
          const dy = (yRange[1] - yRange[0]) / (gridSize - 1);

          const q0 = c * dt;
          const q1 = c2 * dt * dt;
          const q2 = (c * dt / dx) * (c * dt / dx);
          const q3 = (c * dt / dy) * (c * dt / dy);

          const newGrid = new Float32Array(gridSize * gridSize);

          for (let i = 1; i < gridSize - 1; i++) {
            const iOffset = i * gridSize;
            const iPrevOffset = (i - 1) * gridSize;
            const iNextOffset = (i + 1) * gridSize;

            for (let j = 1; j < gridSize - 1; j++) {
              const index = iOffset + j;
              
              if (!this.boundaryMask[index]) {
                newGrid[index] = 0;
                continue;
              }

              const dxx = this.grid[iOffset + j - 1] - 2 * this.grid[index] + this.grid[iOffset + j + 1];
              const dyy = this.grid[iPrevOffset + j] - 2 * this.grid[index] + this.grid[iNextOffset + j];

              const dampingTerm = dampingFactor * q0;
              const denominator = 1.0 + dampingTerm;
              
              newGrid[index] = (
                -this.previousGrid[index] + 
                2 * this.grid[index] + 
                dampingTerm * this.previousGrid[index] +
                q1 * this.sourceGrid[index] +
                q2 * dxx +
                q3 * dyy
              ) / denominator;
            }
          }

          this.previousGrid.set(this.grid);
          this.grid.set(newGrid);
          this.applyBoundaryConditions();
        }

        applyBoundaryConditions() {
          const { gridSize } = this.config;
          for (let i = 0; i < gridSize * gridSize; i++) {
            if (!this.boundaryMask[i]) {
              this.grid[i] = 0;
              this.previousGrid[i] = 0;
            }
          }
        }

        addSourceExcitation(time) {
          const { gridSize } = this.config;
          this.sourceGrid.fill(0);

          if (this.focusI >= 0 && this.focusI < gridSize && this.focusJ >= 0 && this.focusJ < gridSize) {
            const index = this.focusI * gridSize + this.focusJ;
            
            if (this.boundaryMask[index]) {
              // Single pulse at initial time - Gaussian envelope with short duration
              const pulseWidth = 1.0 / this.waveParams.frequency; // One period width
              const pulseDuration = 2.0 * pulseWidth; // Duration of the pulse
              
              let amplitude = 0.0;
              if (time <= pulseDuration) {
                // Gaussian envelope for smooth pulse
                const gaussianWidth = pulseWidth / 3.0;
                const envelope = Math.exp(-Math.pow(time - pulseWidth, 2) / (2 * gaussianWidth * gaussianWidth));
                
                // Single frequency pulse
                amplitude = this.waveParams.amplitude * 10.0 * envelope * 
                           Math.sin(2 * Math.PI * this.waveParams.frequency * time);
              }
              
              this.sourceGrid[index] = amplitude;

              const neighbors = [[-1, 0], [1, 0], [0, -1], [0, 1]];
              for (const [di, dj] of neighbors) {
                const ni = this.focusI + di;
                const nj = this.focusJ + dj;
                
                if (ni >= 0 && ni < gridSize && nj >= 0 && nj < gridSize) {
                  const neighborIndex = ni * gridSize + nj;
                  if (this.boundaryMask[neighborIndex]) {
                    this.sourceGrid[neighborIndex] = amplitude * 0.5;
                  }
                }
              }
            }
          }
        }

        update(dt) {
          this.time += dt;
          this.addSourceExcitation(this.time);
          this.applyWaveEquation(dt);
        }

        getGrid() { return this.grid; }
        getBoundaryMask() { return this.boundaryMask; }
        getTime() { return this.time; }
        reset() {
          this.time = 0;
          this.grid.fill(0);
          this.previousGrid.fill(0);
          this.sourceGrid.fill(0);
        }
      }

      let waveField = null;

      self.addEventListener('message', (event) => {
        const { type, data } = event.data;

        try {
          switch (type) {
            case 'init':
              const { config, waveParams, focusPoint } = data;
              waveField = new WorkerWaveField(config, waveParams, focusPoint);
              
              self.postMessage({
                type: 'initialized',
                data: {
                  grid: waveField.getGrid(),
                  boundaryMask: waveField.getBoundaryMask(),
                  time: waveField.getTime()
                }
              });
              break;

            case 'update':
              if (waveField) {
                const { dt, steps = 1 } = data;
                
                for (let i = 0; i < steps; i++) {
                  waveField.update(dt);
                }

                self.postMessage({
                  type: 'updated',
                  data: {
                    grid: waveField.getGrid(),
                    time: waveField.getTime()
                  }
                });
              }
              break;

            case 'reset':
              if (waveField) {
                waveField.reset();
                
                self.postMessage({
                  type: 'reset',
                  data: {
                    grid: waveField.getGrid(),
                    time: waveField.getTime()
                  }
                });
              }
              break;
          }
        } catch (error) {
          self.postMessage({
            type: 'error',
            data: { message: error.message }
          });
        }
      });
    `;
        return new Blob([workerCode], { type: 'application/javascript' });
    }
    onUpdate(callback) {
        this.updateCallbacks.push(callback);
    }
    async update(dt, steps = 1) {
        if (this.worker && this.isInitialized) {
            this.worker.postMessage({
                type: 'update',
                data: { dt, steps }
            });
        }
    }
    async reset() {
        if (this.worker) {
            this.worker.postMessage({ type: 'reset' });
        }
    }
    getGrid() {
        return this.grid;
    }
    getBoundaryMask() {
        return this.boundaryMask;
    }
    getTime() {
        return this.time;
    }
    isReady() {
        return this.isInitialized;
    }
    destroy() {
        if (this.worker) {
            this.worker.terminate();
            this.worker = null;
        }
        this.updateCallbacks = [];
    }
}
exports.OptimizedWaveField = OptimizedWaveField;
//# sourceMappingURL=optimizedWaveField.js.map