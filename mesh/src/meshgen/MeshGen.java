package meshgen;

import java.io.IOException;
import java.util.ArrayList;

import math.Vector2;
import math.Vector3;

/**
 * A command-line tool for generating, analyzing, and processing meshes stored in the OBJ format.
 * Its three main features are:
 * - Generation of meshes for several basic geometric shapes
 * - Testing meshes for validity and comparing meshes
 * - Generating vertex normals for existing mesh geometry
 *
 * @author ers273, xd93, srm
 * @author Your Name Here
 */

class MeshGen {

    /**
     * Generate the mesh for a cylinder of radius 1 and height 2, centered at the origin.  The
     * cylinder axis is the y axis.  The mesh is generated with <tt>divisions</tt> edges around
     * the top and bottom rims of the cylinder.  The mesh has vertex normals to make the side
     * of the cylinder smooth and the ends flat.  The (u,v) coordinates map the bottom part of
     * the texture to the side of the cylinder and the top part to the ends.  See the assignment
     * writeup for complete details.
     *
     * @param divisions The number of vertices to use to approximate the circular top and bottom rims
     * @return A newly created OBJMesh containing the generated geometry
     */
    public static OBJMesh cylinder(int divisions) {
        OBJMesh outputMesh = new OBJMesh();

        // Task1: Generate Cylinder (20pt)

        float height = 2.0f;
        float radius = 1.0f;
        int topBottomTextureRotationFactor = 0;
        double step = 2 * Math.PI / divisions;
        float textureStep = 1.0f / (float) divisions;
        if (divisions < 2) {
            throw new IllegalArgumentException("divisions for cylinder should be larger than 2");
        }
        // Top center position
        outputMesh.positions.add((new Vector3(0.0f, height / 2, 0.0f))); // 0
        // Bottom center position
        outputMesh.positions.add((new Vector3(0.0f, -height / 2, 0.0f))); // 1
        // Up and down normals
        outputMesh.normals.add((new Vector3(0, 1, 0))); // 0
        outputMesh.normals.add((new Vector3(0, -1, 0))); // 1

        // Top surface
        double phase = -Math.PI / 2.0;
        for (int i = 0; i < divisions; i++) {
            outputMesh.positions.add((new Vector3(radius * (float) Math.cos(i * step + phase), height / 2, radius * (float) Math.sin(i * step + phase))));
            outputMesh.normals.add((new Vector3((float) Math.cos(i * step + phase), 0, (float) Math.sin(i * step + phase))));
            outputMesh.uvs.add((new Vector2(i * textureStep, 0.0f)));
        }
        outputMesh.uvs.add((new Vector2(1.0f, 0.0f)));
        // Bottom surface
        for (int i = 0; i < divisions; i++) {
            outputMesh.positions.add((new Vector3(radius * (float) Math.cos(i * step + phase), -height / 2, radius * (float) Math.sin(i * step + phase))));
            outputMesh.uvs.add((new Vector2(i * textureStep, 0.5f)));
        }
        outputMesh.uvs.add((new Vector2(1.0f, 0.5f)));

        // Top & Bottom texture
        for (int j = 0; j < divisions; j++) {
            int i = (j + topBottomTextureRotationFactor) < divisions ? (j + topBottomTextureRotationFactor) : (j + topBottomTextureRotationFactor) - divisions;
            outputMesh.uvs.add((new Vector2((float) Math.cos(i * step + phase) * 0.25f + 0.75f, (float) -Math.sin(i * step + phase) * 0.25f + 0.75f)));
        }
        for (int j = 0; j < divisions; j++) {
            int i = (j + topBottomTextureRotationFactor) < divisions ? (j + topBottomTextureRotationFactor) : (j + topBottomTextureRotationFactor) - divisions;
            outputMesh.uvs.add((new Vector2((float) -Math.cos(i * step - phase) * 0.25f + 0.25f, (float) -Math.sin(i * step - phase) * 0.25f + 0.75f)));
        }
        outputMesh.uvs.add((new Vector2(0.75f, 0.75f)));
        outputMesh.uvs.add((new Vector2(0.25f, 0.75f)));

        // Top faces
        for (int i = 2; i < divisions + 2; i++) {
            OBJFace triangle = new OBJFace(3, true, true);
            triangle.setVertex(0, 0, 4 * divisions + 2, 0);
            triangle.setVertex(1, (i + 1) < divisions + 2 ? (i + 1) : 2, (2 * divisions + 2) + (i - 1) < 3 * divisions + 2 ? (2 * divisions + 2) + (i - 1) : 2 * divisions + 2, 0);
            triangle.setVertex(2, i, (2 * divisions + 2) + (i - 2), 0);
            outputMesh.faces.add(triangle);
        }
        // Bottom faces
        for (int i = 2 + divisions; i < 2 * divisions + 2; i++) {
            OBJFace triangle = new OBJFace(3, true, true);
            triangle.setVertex(0, 1, 4 * divisions + 3, 1);
            triangle.setVertex(1, i, (2 * divisions + 2) + (i - 2), 1);
            triangle.setVertex(2, (i + 1) < 2 * divisions + 2 ? (i + 1) : 2 + divisions, (2 * divisions + 2) + (i - 1) < 4 * divisions + 2 ? (2 * divisions + 2) + (i - 1) : 3 * divisions + 2, 1);
            outputMesh.faces.add(triangle);
        }
        // Side faces
        for (int i = 2; i < divisions + 2; i++) {
            var positionIndex1 = (i + divisions + 1) < 2 * divisions + 2 ? (i + divisions + 1) : 2 + divisions;
            var positionIndex2 = (i + 1) < divisions + 2 ? (i + 1) : 2;

            OBJFace triangle1 = new OBJFace(3, true, true);
            triangle1.setVertex(0, i, (2 * divisions + 1) - (i - 2), i);
            triangle1.setVertex(1, positionIndex1, (divisions - 1) - (i - 2), positionIndex2);
            triangle1.setVertex(2, i + divisions, divisions - (i - 2), i);
            outputMesh.faces.add(triangle1);

            OBJFace triangle2 = new OBJFace(3, true, true);
            triangle2.setVertex(0, i, (2 * divisions + 1) - (i - 2), i);
            triangle2.setVertex(1, positionIndex2, (2 * divisions) - (i - 2), positionIndex2);
            triangle2.setVertex(2, positionIndex1, (divisions - 1) - (i - 2), positionIndex2);
            outputMesh.faces.add(triangle2);
        }


        return outputMesh;
    }

    /**
     * Generate the mesh for a sphere of radius 1, centered at the origin.  The sphere is triangulated
     * in a longitude-latitude grid, with <tt>divisionsU</tt> divisions around the equator and
     * <tt>divisionsV</tt> rows of triangles from the south to the north pole.  The mesh has the exact
     * surface normals of the geometric sphere, and the (u,v) coordinates are proportional to
     * longitude and latitude.  See the assignment writeup for complete details.
     *
     * @param divisionsU The number of divisions around the equator
     * @param divisionsV The number of divisions from pole to pole
     * @return A newly created OBJMesh containing the generated geometry
     */
    public static OBJMesh sphere(int divisionsU, int divisionsV) {
        OBJMesh outputMesh = new OBJMesh();

        // Task1: Generate Sphere (20pt)
        // TODO:
        double radius = 1.0;
        double equatorStep = 2 * Math.PI / divisionsU;
        double poleStep = Math.PI / divisionsV;
        //add top vertex as 0
        outputMesh.positions.add((new Vector3(0.0f, 1.0f, 0.0f)));
        outputMesh.normals.add((new Vector3(0.0f, 1.0f, 0.0f)));
        int count = 0;
        // Calculate Vertices (positions, uvs, and normals )
        double phase = -Math.PI / 2.0;
        for (int i = 1; i < divisionsV; i++) {
            double phi = poleStep * i;
            for (int j = 0; j < divisionsU; j++) {
                count++;
                outputMesh.positions.add((new Vector3(
                        (float) (radius * Math.sin(phi) * Math.cos(j * equatorStep + phase)),
                        (float) (radius * Math.cos(phi)),
                        (float) (radius * Math.sin(phi) * Math.sin(j * equatorStep + phase))
                )));
                outputMesh.normals.add((new Vector3(
                        (float) (radius * Math.sin(phi) * Math.cos(j * equatorStep + phase)),
                        (float) (radius * Math.cos(phi)),
                        (float) (radius * Math.sin(phi) * Math.sin(j * equatorStep + phase))
                )));
            }
        }
        outputMesh.positions.add(((new Vector3(0.0f, -1.0f, 0.0f))));
        outputMesh.normals.add((new Vector3(0.0f, -1.0f, 0.0f)));

        // Calculate texture
        int totalTextureNum = (divisionsU + 1) * (divisionsV + 1);
        float equatorTextureStep = 1.0f / (float) divisionsU;
        float poleTextureStep = 1.0f / (float) divisionsV;
        for (int i = 0; i <= divisionsV; i++) {
            float poleTextureCoordinate = i * poleTextureStep;
            for (int j = 0; j <= divisionsU; j++) {
                float equatorTextureCoordinate = j * equatorTextureStep;
                outputMesh.uvs.add((new Vector2(equatorTextureCoordinate, poleTextureCoordinate)));
            }

        }

        // Calculate indices in faces (use OBJFace class)
        for (int i = 0; i < divisionsV - 1; i++) {
            for (int j = 1; j < divisionsU + 1; j++) {
                //draw up triangle
                OBJFace triangle = new OBJFace(3, true, true);
                triangle.setVertex(0, Math.max(0, j == divisionsU ? 1 + (i - 1) * divisionsU : j + 1 + (i - 1) * divisionsU), (totalTextureNum - 1) - i * (divisionsU + 1) - j, Math.max(0, j == divisionsU ? 1 + (i - 1) * divisionsU : j + 1 + (i - 1) * divisionsU));
                triangle.setVertex(1, j == divisionsU ? (1 + i * divisionsU) : j + i * divisionsU + 1, (totalTextureNum - 1) - (i + 1) * (divisionsU + 1) - j, j == divisionsU ? (1 + i * divisionsU) : j + i * divisionsU + 1);
                triangle.setVertex(2, j + i * divisionsU, (totalTextureNum - 1) - (i + 1) * (divisionsU + 1) - (j - 1), j + i * divisionsU);

                outputMesh.faces.add(triangle);
                OBJFace triangleDown = new OBJFace(3, true, true);
                triangleDown.setVertex(0, j + i * divisionsU, (totalTextureNum - 1) - (i + 1) * (divisionsU + 1) - (j - 1), j + i * divisionsU);
                triangleDown.setVertex(1, j == divisionsU ? (1 + i * divisionsU) : j + i * divisionsU + 1, (totalTextureNum - 1) - (i + 1) * (divisionsU + 1) - j, j == divisionsU ? (1 + i * divisionsU) : j + i * divisionsU + 1);
                triangleDown.setVertex(2, Math.min((divisionsV - 1) * divisionsU + 1, j + (i + 1) * divisionsU), (totalTextureNum - 1) - (i + 2) * (divisionsU + 1) - (j - 1), Math.min((divisionsV - 1) * divisionsU + 1, j + (i + 1) * divisionsU));
                outputMesh.faces.add(triangleDown);
            }
        }
        return outputMesh;
    }


    /**
     * Create vertex normals for an existing mesh.  The triangles, positions, and uvs are copied
     * from the input mesh, and new vertex normals are computed by averaging the normals of the
     * faces that share each vertex.
     *
     * @param inputMesh The input mesh, whose triangles and vertex positions define the normals
     * @return A newly created OBJMesh that is a copy of the input mesh but with new normals
     */
    public static OBJMesh createNormals(OBJMesh inputMesh) {
        OBJMesh outputMesh = new OBJMesh();

        // Task2: Compute Normals (35pt)
        // TODO:
        // Copy position data
        for (Vector3 vec : inputMesh.positions) {
            outputMesh.positions.add(vec.clone());
        }
        // Copy UV data
        for (Vector2 vec : inputMesh.uvs) {
            outputMesh.uvs.add(vec.clone());
        }
        // Each vertex gets a unique normal
        for (Vector3 vec : inputMesh.positions) {
            outputMesh.normals.add((new Vector3(0.f, 0.f, 0.f)));
        }
        for (OBJFace face : inputMesh.faces) {
            outputMesh.faces.add(face);
            Vector3 a = outputMesh.positions.get(face.positions[0]).clone();
            Vector3 b = outputMesh.positions.get(face.positions[1]).clone();
            Vector3 c = outputMesh.positions.get(face.positions[2]).clone();
            Vector3 line1 = a.clone().sub(b.clone());
            Vector3 line2 = a.clone().sub(c.clone());
            Vector3 direction = line1.cross(line2).normalize();
            outputMesh.normals.set(face.positions[0], outputMesh.normals.get(face.positions[0]).clone().add(direction.clone()));
            outputMesh.normals.set(face.positions[1], outputMesh.normals.get(face.positions[1]).clone().add(direction.clone()));
            outputMesh.normals.set(face.positions[2], outputMesh.normals.get(face.positions[2]).clone().add(direction.clone()));
        }
        for (int i = 0; i < outputMesh.normals.size(); i++) {
            Vector3 tempNorm = outputMesh.normals.get(i).clone();
            outputMesh.normals.set(i, tempNorm.normalize().clone());
        }
        // Initialize output faces
        // Calculate face normals, distribute to adjacent vertices
        // Normalize new normals
        for (OBJFace face : outputMesh.faces) {
            face.normals = face.positions.clone();
        }

        return outputMesh;
    }


    //
    // The following are extra credits, it is not required, do as you are interested in it.
    //

    /**
     * Generate the mesh for a torus of major radius 1 and minor radius <tt>minorRadius</tt>.
     * The symmetry axis is the y axis.  The torus is triangulated in a grid with
     * <tt>divisionsU</tt> divisions along the tube (major direction) and <tt>divisionsV</tt>
     * divisions around the tube (minor direction).  The vertex normals are the exact normals
     * to the geometric torus, and the (u,v) coordinates follow the triangulation grid.
     * See the assignment writeup for complete details.
     *
     * @param divisionsU  The number of divisions in the major direction
     * @param divisionsV  The number of divisions in the minor direction
     * @param minorRadius The minor radius (radius of the tube)
     * @return A newly created OBJMesh containing the generated geometry
     */
    public static OBJMesh torus(int divisionsU, int divisionsV, float minorRadius) {
        OBJMesh outputMesh = new OBJMesh();

        // Extra Credit: Generate Turos (10pt)
        // TODO:
        float majorRadius = 1.0f;
        double thetaStep = 2 * Math.PI / divisionsV;
        double phiStep = 2 * Math.PI / divisionsU;
        double phase = -Math.PI / 2.0;
        double textureStepU = 1.0f / (float) divisionsU;
        double textureStepV = 1.0f / (float) divisionsV;

        // Calculate vertices: positions, uvs and normals
        for (int i = 0; i < divisionsU; i++) {
            for (int j = 0; j < divisionsV; j++) {
                float xCoordinate = (float) ((majorRadius - minorRadius * Math.cos(j * thetaStep)) * (Math.cos(i * phiStep + phase)));
                float zCoordinate = (float) ((majorRadius - minorRadius * Math.cos(j * thetaStep)) * (Math.sin(i * phiStep + phase)));
                float yCoordinate = (float) (minorRadius * Math.sin(j * thetaStep));
                outputMesh.positions.add((new Vector3(xCoordinate, yCoordinate, zCoordinate)));
                outputMesh.normals.add(new Vector3((float) (-minorRadius * Math.cos(j * thetaStep) * Math.cos(i * phiStep + phase)),
                        (float) (minorRadius * Math.sin(j * thetaStep)),
                        (float) (-minorRadius * Math.cos(j * thetaStep) * Math.sin(i * phiStep + phase))));

            }
        }

        for (int i = 0; i < divisionsU + 1; i++) {
            for (int j = 0; j < divisionsV + 1; j++) {
                    outputMesh.uvs.add(new Vector2((float)(1-i*textureStepU),(float)(1 - j*textureStepV)));
            }
        }

        // Calculate indices on faces (use OBJFace class)
        for (int i = 0; i < divisionsU; i++) {
            for (int j = 0; j < divisionsV; j++) {
                int jj = j + 1;
                int ii = i - 1;
                int itop = i + 1;
                if (j == divisionsV - 1) {
                    jj = 0;
                }
                if (i == 0) {
                    ii = divisionsU - 1;
                }
                if (i == divisionsU - 1) {
                    itop = 0;
                }
                int sssi = i;
                if (i == 0) {
                    sssi = divisionsU ;
                }
                OBJFace topTriangle = new OBJFace(3, true, true);
                topTriangle.setVertex(0, j + i * divisionsV, j + sssi * (divisionsV + 1), j + i * divisionsV);
                topTriangle.setVertex(2, jj + ii * divisionsV, j + 1 + (i == 0 ? (divisionsU - 1) * (divisionsV + 1) : (i-1) * (divisionsV + 1)), jj + ii * divisionsV);
                topTriangle.setVertex(1, jj + i * divisionsV, j + 1 + sssi * (divisionsV + 1), jj + i * divisionsV);
                outputMesh.faces.add(topTriangle);
                if (i == 0) {
                    System.out.println("v1: " + String.valueOf(j + i * (divisionsV + 1)) + " v2: " + String.valueOf(j + 1 + (divisionsU - 1) * (divisionsV + 1) ) + " v 3: " + String.valueOf(j + 1 + i * (divisionsV + 1)));
                }

                OBJFace bottomTriangle = new OBJFace(3, true, true);
                bottomTriangle.setVertex(0, j + i * divisionsV, j + i * (divisionsV + 1), j + i * divisionsV);
                bottomTriangle.setVertex(1, j + itop * divisionsV, j + (i + 1) * (divisionsV + 1), j + itop * divisionsV);
                bottomTriangle.setVertex(2, jj + i * divisionsV, j + 1 + i * (divisionsV + 1), jj + i * divisionsV);
                outputMesh.faces.add(bottomTriangle);

            }
        }


        return outputMesh;
    }


    public static OBJMesh geodesicSphere(int divisionU, int divisionV) {
        OBJMesh outputMesh = new OBJMesh();

        // Extra Credit: Geodesic Sphere (10pt)
        // TODO:
        // Calculate vertices: positions, uvs and normals
        // Calculate indices on faces (use OBJFace class)

        return outputMesh;
    }

    public static void main(String[] args) {
        if (args == null || args.length < 2) {
            System.err.println("Error: not enough input arguments.");
            printUsage();
            System.exit(1);
        }

        if (args[0].equals("-g")) { // Generate mesh
            int divisionsU = 32;
            int divisionsV = 16;
            float minorRadius = 0.25f;
            String outputFilename = null;

            for (int i = 2; i < args.length; i += 2) {
                if (i + 1 == args.length) {
                    System.err.println("Error: expected argument after \"" + args[i] + "\" flag.");
                    printUsage();
                    System.exit(1);
                }
                if (args[i].equals("-n")) { // Divisions latitude
                    divisionsU = Integer.parseInt(args[i + 1]);
                } else if (args[i].equals("-m")) { // Divisions longitude
                    divisionsV = Integer.parseInt(args[i + 1]);
                } else if (args[i].equals("-r")) { // Inner radius
                    minorRadius = Float.parseFloat(args[i + 1]);
                } else if (args[i].equals("-o")) { // Output filename
                    outputFilename = args[i + 1];
                } else {
                    System.err.println("Error: Unknown option \"" + args[i] + "\"");
                    printUsage();
                    System.exit(1);
                }
            }

            if (outputFilename == null) {
                System.err.println("Error: expected -o argument.");
                printUsage();
                System.exit(1);
            }

            OBJMesh outputMesh = null;
            if (args[1].equals("cylinder")) {
                outputMesh = cylinder(divisionsU);
            } else if (args[1].equals("sphere")) {
                outputMesh = sphere(divisionsU, divisionsV);
            } else if (args[1].equals("torus")) {
                outputMesh = torus(divisionsU, divisionsV, minorRadius);
            } else {
                System.err.println("Error: expected geometry type.");
                printUsage();
                System.exit(1);
            }

            System.out.println("Output mesh is valid: " + outputMesh.isValid(true));

            try {
                outputMesh.writeOBJ(outputFilename);
            } catch (IOException e) {
                System.err.println("Error: could not write file " + outputFilename);
                System.exit(1);
            }

        } else if (args[0].equals("-i")) { // Assign normals
            if (!args[2].equals("-o")) {
                System.err.println("Error: expected -o argument.");
                printUsage();
                System.exit(1);
            }
            OBJMesh inputMesh = null;
            try {
                inputMesh = new OBJMesh(args[1]);
            } catch (OBJMesh.OBJFileFormatException e) {
                System.err.println("Error: Malformed input OBJ file: " + args[1]);
                System.err.println(e);
                System.exit(1);
            } catch (IOException e) {
                System.err.println("Error: could not read file " + args[1]);
                System.err.println(e);
                System.exit(1);
            }
            OBJMesh outputMesh = createNormals(inputMesh);

            System.out.println("Output mesh is valid: " + outputMesh.isValid(true));

            try {
                outputMesh.writeOBJ(args[3]);
            } catch (IOException e) {
                System.err.println("Error: could not write file " + args[3]);
                System.exit(1);
            }

        } else if (args[0].equals("-v")) { // Verify an OBJ file
            if (args.length != 2) {
                System.err.println("Error: expected an input file argument.");
                printUsage();
                System.exit(1);
            }
            OBJMesh mesh = null;
            try {
                mesh = new OBJMesh(args[1]);
            } catch (OBJMesh.OBJFileFormatException e) {
                System.err.println("Error: Malformed input OBJ file: " + args[1]);
                System.err.println(e);
                System.exit(1);
            } catch (IOException e) {
                System.err.println("Error: could not read file " + args[1]);
                System.err.println(e);
                System.exit(1);
            }
            System.out.println("Input mesh is valid OBJ syntax: " + mesh.isValid(true));

        } else if (args[0].equals("-c")) { // Compare two OBJ files
            float eps = 1e-5f;
            if (args.length != 3 && args.length != 5) {
                System.err.println("Error: expected 2 input file arguments and optional epsilon.");
                printUsage();
                System.exit(1);
            }
            if (args.length == 5) {
                if (!args[1].equals("-e")) {
                    System.err.println("Error: expected -e flag after -c.");
                    printUsage();
                    System.exit(1);
                }
                eps = Float.parseFloat(args[2]);
            }
            OBJMesh m1 = null, m2 = null;
            int m1arg = (args.length == 3) ? 1 : 3;
            int m2arg = (args.length == 3) ? 2 : 4;
            try {
                m1 = new OBJMesh(args[m1arg]);
            } catch (OBJMesh.OBJFileFormatException e) {
                System.err.println("Error: Malformed input OBJ file: " + args[m1arg]);
                System.err.println(e);
                System.exit(1);
            } catch (IOException e) {
                System.err.println("Error: could not read file " + args[m1arg]);
                System.err.println(e);
                System.exit(1);
            }
            try {
                m2 = new OBJMesh(args[m2arg]);
            } catch (OBJMesh.OBJFileFormatException e) {
                System.err.println("Error: Malformed input OBJ file: " + args[m2arg]);
                System.err.println(e);
                System.exit(1);
            } catch (IOException e) {
                System.err.println("Error: could not read file " + args[m2arg]);
                System.err.println(e);
                System.exit(1);
            }
            System.out.println("Meshes are equivalent: " + OBJMesh.compare(m1, m2, true, eps));

        } else {
            System.err.println("Error: Unknown option \"" + args[0] + "\".");
            printUsage();
            System.exit(1);
        }
    }

    public static void printUsage() {
        System.out.println("Usage:");
        System.out.println("(1) MeshGen -g <cylinder|sphere|torus> [-n <divisionsLatitude>] [-m <divisionsLongitude>] [-r <innerRadius>] -o <output.obj>");
        System.out.println("(2) MeshGen -i <input.obj> -o <output.obj>");
        System.out.println("(3) MeshGen -v <input.obj>");
        System.out.println("(4) MeshGen -c [-e <epsilon>] <m1.obj> <m2.obj>");
        System.out.println();
        System.out.println("(1) creates an OBJ mesh of a cylinder, sphere, or torus.");
        System.out.println("Cylinder ignores -n and -r flags, sphere ignores the -r flag.");
        System.out.println("(2) takes in an OBJ mesh, strips it of normals (if any), and assigns new normals based on the normals of its faces.");
        System.out.println("(3) verifies that an input OBJ mesh file conforms to the OBJ standard.");
        System.out.println("(4) compares two input OBJ files and checks if they are equivalent up to an optional epsilon parameter (by default, epsilon=1e-5).");
    }

}
