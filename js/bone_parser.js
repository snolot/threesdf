const BoneParser = (raymarcher, modelURL) => {
	
	const _scene = new THREE.Scene();
	const _loader = new THREE.GLTFLoader();
	const _clock = new THREE.Clock();
	const _rm = raymarcher;
    
	const _bonesNames = [
		'mixamorigHead',
		'mixamorigNeck',
		'mixamorigLeftShoulder',
		'mixamorigLeftForeArm',
		'mixamorigLeftHand',
		'mixamorigRightShoulder',
		'mixamorigRightForeArm',
		'mixamorigRightHand',
		'mixamorigSpine',
		'mixamorigLeftUpLeg',
		'mixamorigLeftLeg',
		'mixamorigLeftFoot',
		'mixamorigRightUpLeg',
		'mixamorigRightLeg',
		'mixamorigRightFoot',
	];
	let _bonesRef = [];
    let _paused = false;
    let _markers = [];
    let _vertices = [];
    //let _values = [];
    let _speed = 1;
    let _walkertime = 0;
    let _initphase = 0;
    let _sizes = [];
    let _lastTime = Date.now() * .001;
    let _startTime = _lastTime;
    let _mixer;
    let _ready = false;
    let _notCalibratedOnTarget = false;
    let _material;
    
    for ( i = 0; i < 15; i++ ){
        _vertices.push( new THREE.Vector3(0,0,0) );
    	//_values.push( 0,0,0 );
    }

    _rm.setUniform( "anchors", "v3v", _vertices );

    const _reset = _ => {
    	_speed	= 1;
        _walkertime = 0;
        _initphase = 0;
    }

    _loader.load(modelURL, gltf => {
		console.log('model loaded.');
		const model = gltf.scene;
		const animations = gltf.animations;
		_rm.addToScene(model);

		_mixer = new THREE.AnimationMixer( model );
		var action = _mixer.clipAction( animations[ 0 ] );
		action.play();
		
		_bonesNames.map(name => {
			console.log(name)
			_bonesRef.push(model.getObjectByName(name));
		})

		_material = _rm.getMaterial();
		_ready = true;
	});

    return {
    	update: time => {
    		if(!_ready) return;

    		const delta = _clock.getDelta();

			if(_mixer){
				_mixer.update( delta );
			}

    		const curtime = Date.now() * .0005 - _startTime;
    		const timediff = (curtime - _lastTime);

    		i = 0;
        	//_values = [];
        	_sizes = [];
	        
	        const bonePosition = new THREE.Vector3();
	        const _target = _rm.getTarget();

	        while ( i < 15 ){
	        	const bone = _bonesRef[i];
	        	
	        	bonePosition.setFromMatrixPosition(bone.matrixWorld);//.add(bone.position);

	        	_vertices[ i ].copy(bonePosition);
	            
	        	var v = _vertices[i].multiplyScalar( 10 );
	            //_values.push( v.x, v.y, v.z);

	            if( i == 8/* && _notCalibratedOnTarget*/){
	            	_notCalibratedOnTarget = true;
	                _target.copy( v );
	                _target.y += 2;
	            }

	            i++;

	        }

        	_material.uniforms.anchors.value = _vertices;
        	_material.uniforms.anchors.needsUpdate = true;

        	_lastTime = curtime;
    	}
    }
}