"use strict";
Object.defineProperty(exports, "__esModule", { value: true });
exports.Parabola = void 0;
class Parabola {
    constructor(diameter, focus, vertex, concaveUp = true) {
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
    getY(x) {
        const { coefficient, vertex } = this.params;
        return coefficient * Math.pow(x - vertex.x, 2) + vertex.y;
    }
    /**
     * Check if a point is inside the parabolic boundary
     */
    isPointInside(point) {
        const { diameter, vertex, concaveUp } = this.params;
        // Check if within horizontal bounds
        if (Math.abs(point.x - vertex.x) > diameter / 2) {
            return false;
        }
        const parabolaY = this.getY(point.x);
        if (concaveUp) {
            return point.y >= parabolaY;
        }
        else {
            return point.y <= parabolaY;
        }
    }
    /**
     * Get the normal vector at a given x position on the parabola
     */
    getNormalVector(x) {
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
    reflect(incidentPoint, incidentDirection) {
        const normal = this.getNormalVector(incidentPoint.x);
        // Reflection formula: r = d - 2(d·n)n
        const dotProduct = incidentDirection.x * normal.x + incidentDirection.y * normal.y;
        return {
            x: incidentDirection.x - 2 * dotProduct * normal.x,
            y: incidentDirection.y - 2 * dotProduct * normal.y
        };
    }
    getParams() {
        return { ...this.params };
    }
}
exports.Parabola = Parabola;
//# sourceMappingURL=parabola.js.map