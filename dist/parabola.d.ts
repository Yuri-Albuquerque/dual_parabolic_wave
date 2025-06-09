import { Point2D, ParabolaParams } from './types';
export declare class Parabola {
    private params;
    constructor(diameter: number, focus: number, vertex: Point2D, concaveUp?: boolean);
    /**
     * Calculate y-coordinate for given x on the parabola
     * Equation: y = a(x - h)² + k, where (h,k) is vertex
     */
    getY(x: number): number;
    /**
     * Check if a point is inside the parabolic boundary
     */
    isPointInside(point: Point2D): boolean;
    /**
     * Get the normal vector at a given x position on the parabola
     */
    getNormalVector(x: number): Point2D;
    /**
     * Calculate reflection point for a ray hitting the parabola
     */
    reflect(incidentPoint: Point2D, incidentDirection: Point2D): Point2D;
    getParams(): ParabolaParams;
}
//# sourceMappingURL=parabola.d.ts.map