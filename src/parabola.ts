import { Point2D, ParabolaParams } from './types';

export class Parabola {
  private params: ParabolaParams;

  constructor(diameter: number, focus: number, vertex: Point2D, concaveUp: boolean = true) {
    this.params = {
      diameter,
      focus,
      vertex,
      coefficient: concaveUp ? 1 / (4 * focus) : -1 / (4 * focus),
      concaveUp
    };
  }

  /**
   * Calculate y-coordinate for given x on the parabola
   * Equation: y = a(x - h)² + k, where (h,k) is vertex
   */
  getY(x: number): number {
    const { coefficient, vertex } = this.params;
    return coefficient * Math.pow(x - vertex.x, 2) + vertex.y;
  }

  /**
   * Check if a point is inside the parabolic boundary
   */
  isPointInside(point: Point2D): boolean {
    const { diameter, vertex, concaveUp } = this.params;
    
    // Check if within horizontal bounds
    if (Math.abs(point.x - vertex.x) > diameter / 2) {
      return false;
    }

    const parabolaY = this.getY(point.x);
    
    if (concaveUp) {
      return point.y >= parabolaY;
    } else {
      return point.y <= parabolaY;
    }
  }

  /**
   * Get the normal vector at a given x position on the parabola
   */
  getNormalVector(x: number): Point2D {
    const { coefficient, vertex } = this.params;
    
    // Derivative: dy/dx = 2a(x - h)
    const slope = 2 * coefficient * (x - vertex.x);
    
    // Normal vector is perpendicular to tangent
    const normalMag = Math.sqrt(1 + slope * slope);
    return {
      x: -slope / normalMag,
      y: 1 / normalMag
    };
  }

  /**
   * Calculate reflection point for a ray hitting the parabola
   */
  reflect(incidentPoint: Point2D, incidentDirection: Point2D): Point2D {
    const normal = this.getNormalVector(incidentPoint.x);
    
    // Reflection formula: r = d - 2(d·n)n
    const dotProduct = incidentDirection.x * normal.x + incidentDirection.y * normal.y;
    
    return {
      x: incidentDirection.x - 2 * dotProduct * normal.x,
      y: incidentDirection.y - 2 * dotProduct * normal.y
    };
  }

  getParams(): ParabolaParams {
    return { ...this.params };
  }
}
