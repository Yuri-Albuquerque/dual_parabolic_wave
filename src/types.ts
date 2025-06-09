export interface Point2D {
  x: number;
  y: number;
}

export interface ParabolaParams {
  diameter: number;
  focus: number;
  vertex: Point2D;
  coefficient: number;
  concaveUp: boolean;
}

export interface WaveParams {
  frequency: number;
  wavelength: number;
  speed: number;
  amplitude: number;
}

export interface SimulationConfig {
  gridSize: number;
  xRange: [number, number];
  yRange: [number, number];
  timeStep: number;
  dampingFactor: number;
  reflectionCoeff: number;
}
