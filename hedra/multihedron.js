import * as THREE from 'three.min.js';
import { OrbitControls } from 'OrbitControls.js';
import { ConvexGeometry } from 'ConvexGeometry.js';
import { GLTFExporter } from 'GLTFExporter.js';

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, 1000 );  //default 75 fov
const controls = new OrbitControls( camera, document.getElementById("bg") );
camera.position.z = 25;

const renderer = new THREE.WebGLRenderer({canvas: document.getElementById("bg"), alpha: true});
renderer.setSize( window.innerWidth, window.innerHeight );
renderer.setPixelRatio( window.devicePixelRatio );

const stdMaterial2 = new THREE.MeshStandardMaterial( { color: 0x00FF00} ); //reacts to light

//Multihedron
const hedronPoints = [{x:1,y:2,z:3},{x:0.5,y:2,z:3},{x:0.5,y:1,z:3},{x:0.5,y:1,z:1.5},{x:2,y:1,z:3},{x:2,y:0.5,z:3},{x:1,y:0.5,z:3},{x:1,y:0.5,z:1.5},{x:3,y:1,z:2},{x:3,y:0.5,z:2},{x:3,y:0.5,z:1},{x:1.5,y:0.5,z:1},{x:3,y:2,z:1},{x:3,y:2,z:0.5},{x:3,y:1,z:0.5},{x:1.5,y:1,z:0.5},{x:1,y:4,z:1},{x:0.5,y:4,z:1},{x:1,y:4,z:0.5},{x:0.5,y:4,z:0.5},{x:0.5,y:2,z:0.5}];
for(var i = 0; i<hedronPoints.length; i++){
  hedronPoints[i] = new THREE.Vector3(hedronPoints[i].x, hedronPoints[i].y, hedronPoints[i].z);
}
const multihedronGeo = new ConvexGeometry( hedronPoints );
const multihedron = new THREE.Mesh( multihedronGeo, stdMaterial2 );
multihedron.position.z = 0;
multihedron.position.y = 0;
scene.add( multihedron );

//var lineGeometry = multihedronGeo.clone();
//  lineGeometry.scale(1.0001, 1.0001, 1.0001);
//const edges = new THREE.EdgesGeometry( lineGeometry );
const edges = new THREE.EdgesGeometry( multihedronGeo );
const line = new THREE.LineSegments( edges, new THREE.LineBasicMaterial( { color: 0xffffff } ) );
  //line.position.addScalar(-0.0001, -0.0001, -0.0001);
scene.add( line );

const hemiLight = new THREE.HemisphereLight( 0xffffff, 0xffffff, 0.6 );
hemiLight.color.setHSL( 0.6, 1, 0.6 );
hemiLight.groundColor.setHSL( 0.095, 1, 0.75 );
hemiLight.position.set( 0, 50, 0 );
scene.add( hemiLight );
/*
const ambientLight = new THREE.AmbientLight(0x555555);
scene.add( ambientLight );

const pointLight = new THREE.PointLight(0xAAAAAA);
pointLight.position.set(0,10,0);
scene.add( pointLight );

//Shows Position of point light
const lightHelper = new THREE.PointLightHelper(pointLight);
scene.add(lightHelper);
*/
function animate(){
  renderer.render( scene, camera );
  requestAnimationFrame(animate);
}

// Instantiate a exporter
const exporter = new GLTFExporter();

// Parse the input and generate the glTF output
exporter.parse(
	multihedron,
	// called when the gltf has been generated
	function ( gltf ) {

    //copy output to a .gltf file
		console.log( JSON.stringify(gltf) );
    //var file = new File(gltf.buffers, "associahedron.glb");
    //console.log(file)
    //window.location.href = gltf.buffers[0].uri;
    //var element = document.createElement('a');
    //element.setAttribute('href', 'gltf.buffers[0].uri');
    //element.setAttribute('download', "multihedron.glb");

    //element.style.display = 'none';
    //document.body.appendChild(element);

	},
	// called when there is an error in the generation
	function ( error ) {

		console.log( 'An error happened' );

	},
	//{"binary": true}
);




requestAnimationFrame(animate);

window.addEventListener("resize", function(){
  renderer.setSize( window.innerWidth, window.innerHeight );
  renderer.setPixelRatio( window.devicePixelRatio );
  camera.aspect = window.innerWidth / window.innerHeight;
  camera.updateProjectionMatrix();
});
