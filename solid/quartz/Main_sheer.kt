import IOData.IOGroXtc.readStructureFromGRO
import IOData.IOGroXtc.writeStructureToGRO
import IOData.IOItp.readMoleculeFromITP
import IOData.packageParameters.Companion.maximumDistance
import IOData.packageParameters.Companion.minimalDistance
import IOData.packageParameters.Companion.packageParameter
import IOData.packageParameters.Companion.readInputParametersFromFile
import IOData.systemTopology.molPathGROList
import IOData.systemTopology.molPathITPList
import IOData.systemTopology.molTypesNumber
import IOData.systemTopology.readMoleculesTypesAndOrder
import IOData.systemTopology.systemITP
import assembler.Assembler.findPositionByRotation
import assembler.Assembler.insertCOMIntoShape
import assembler.Assembler.wedgeOffParallel
import geometry.*
import methods.progressBar.drawProgressBar
import java.lang.Math.*
import java.util.*
import kotlin.system.exitProcess

fun main(args: Array<String>) {
  Locale.setDefault(Locale.US)
  readInputParametersFromFile("package_parameters.ass")
  readMoleculesTypesAndOrder("systems.list")
  // BASE STRUCTURE
  val adsorbent = "qua" //  "car" //""sil" // "sil" // "gr"
  var slitWidth = 6.0 //nm
  val delta = 0.6 // nm
  slitWidth -= 2*delta

  val octane = intArrayOf(457)
  val water =  intArrayOf(12119)



  val lx = 7.970756
  val ly = 10.354234
  val lz = 9.162737

  val octBorders = Point3D(3.991534, 10.35423, 3.004634)

  val size = Point3D(lx, ly, lz)
  val solid = size.z - slitWidth + delta

  val maxSheer = -5.978207
  val sheer = 1.0/sqrt(3.0)

  for (i in octane.indices) { // number of systems
    val walls = readStructureFromGRO("${adsorbent}.gro")
    for (atom in walls.atomsAsArrayList) {
      atom.moleculeNumber = 0
    }
    var structure = walls.copy()
    structure.size = size
    structure = structure.applyStructurePBC()
    writeStructureToGRO("qua_cubic.gro", structure)
    exitProcess(0)
    /**
    structure.atomsAsArrayList.indices.forEach {
      val at = structure.getAtom(it)
      val localSheer = at.coordinates.y*sheer
      at.coordinates = at.coordinates.plus(Point3D(localSheer, 0.0 , 0.0))
    }
    */

    val newCentersList = structure.allAtomCoordinates//new PointList();
    val molNamesList = ArrayList<String>()
    var atomNumber = walls.numberOfAtoms

    val center = structure.center
    //val droplet = Sphere(center, 2.0)
    //val sneaky = Parallelepiped(center, Point3D(3.0, 3.0, 3.0)) // only for small C of phobic liquid
    val oct = Parallelepiped(Point3D(center.x, center.y, solid + slitWidth/4),
        octBorders.minus(Point3D(0.0,0.0,delta)))
    val wat1 = Parallelepiped(Point3D((size.x - octBorders.x)/4, center.y, solid + slitWidth/2),
        Point3D((size.x - octBorders.x)/2, size.y, slitWidth))
    val wat2 = Parallelepiped(Point3D((size.x - octBorders.x)/4*3 + octBorders.x, center.y, solid + slitWidth/2),
        Point3D((size.x - octBorders.x)/2, size.y, slitWidth))
    //val wat3 = Parallelepiped(Point3D(center.x, (size.y - octBorders.y)/4, solid + slitWidth/2),
     //   Point3D(octBorders.x, (size.y - octBorders.y)/2, slitWidth))
    //val wat4 = Parallelepiped(Point3D(center.x, (size.y - octBorders.y)/4*3 + octBorders.y, solid + slitWidth/2),
    //    Point3D(octBorders.x, (size.y - octBorders.y)/2, slitWidth))
    val wat5 = Parallelepiped(Point3D(center.x, center.y, solid + 3.0*slitWidth/4.0), octBorders)
    val shapes = listOf(oct, wat1, wat2, wat5)
    val totalVolume = wat1.getVolume() +
        wat2.getVolume() +
        //wat3.getVolume() +
        //wat4.getVolume() +
        wat5.getVolume()
    val w1 = round(water[i]*wat1.getVolume()/totalVolume).toInt()
    val w2 = round(water[i]*wat2.getVolume()/totalVolume).toInt()
    //val w3 = round(water[i]*wat3.getVolume()/totalVolume).toInt()
    //val w4 = round(water[i]*wat4.getVolume()/totalVolume).toInt()
    var w5 = round(water[i]*wat5.getVolume()/totalVolume).toInt()
    val totalWat = w1+ w2 + w5
    w5 = w5 + water[i] - totalWat

    val molNumbersList = arrayListOf(
        octane[i],
        w1,
        w2,
      //  w3,
      //  w4,
        w5
    )
    //val outFileName = "${adsorbent}_${molName}_dw_${surfactant[i]}_${decane[i]}_${water[i]}.gro"
    val outFileName = "${adsorbent}_ow_${octane[i]}_${water[i]}.gro"

    for (k in 0 until molTypesNumber) {
      val currentMolecule = readStructureFromGRO(molPathGROList[k]) // gro
      val currentTopology = readMoleculeFromITP(molPathITPList[k]) // top
      molNamesList.add(currentTopology.name)
      println("Insertion of " + molNamesList[k] + " in progress")
      systemITP.append(molNamesList[k]).append(" ").append(molNumbersList[k]).append("\n")
      val oldAtomsList = currentMolecule.allAtomCoordinates.centerPoints()
      for (j in 0 until molNumbersList[k]) { // insert molecules
        drawProgressBar(j, molNumbersList[k])
        val newCenterOfMass = insertCOMIntoShape(
            j,
            structure.size,
            oldAtomsList.moleculeSize(),
            newCentersList, 
            shapes[k],
            packageParameter
        )
        newCentersList.addPoint(newCenterOfMass)
        currentMolecule.allAtomCoordinates = findPositionByRotation(
            minimalDistance,
            structure,
            oldAtomsList,
            newCenterOfMass,
            1e2.toInt())
        for (localNumber in 0 until currentMolecule.numberOfAtoms) {
          val atom = currentMolecule.getAtom(localNumber).copy()
          atom.number = atomNumber
          atom.moleculeNumber = newCentersList.numberOfPoints - walls.numberOfAtoms
          structure.addAtom(atom)
          atomNumber++
        }
      }
    }
    println((structure.getAtom(structure.numberOfAtoms - 1).moleculeNumber + 1).toString() +
        " molecules in the box")
    structure = structure.applyStructurePBC()
    structure = wedgeOffParallel(minimalDistance, maximumDistance, structure)

    structure.atomsAsArrayList.indices.forEach {
      val at = structure.getAtom(it)
      val localSheer = at.coordinates.y*sheer
      at.coordinates = at.coordinates.plus(Point3D(-localSheer, 0.0 , 0.0))
    }

    writeStructureToGRO(outFileName, structure)
    println("Structure $outFileName ready")

    //val currentSystemITP = BufferedWriter(FileWriter(File("system.itp")))
    //currentSystemITP.write(systemITP.toString())
    //currentSystemITP.close()
    //println("System.itp also ready")
    //println()
  }
}
