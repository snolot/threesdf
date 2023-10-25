const RayMarcher = _ => {
	
	const _tl = new THREE.TextureLoader();
    const _cl = new THREE.CubeTextureLoader();
    const _rc = new THREE.Raycaster();
    const _mouse = new THREE.Vector2();
    const _scene = new THREE.Scene();
    const _renderer = new THREE.WebGLRenderer({powerPreference: "high-performance", transparent:false});
    const _renderCamera = new THREE.OrthographicCamera(-1,1,1,-1,1/Math.pow( 2, 53 ),1);
    const _geom = new THREE.BufferGeometry();
    const _mesh = new THREE.Mesh( _geom, null );
    const _camera = new THREE.PerspectiveCamera( 60, 1, 0.1,1 );
    const _target = new THREE.Vector3();

	let _distance;
	let _precision;
	let _loaded = false;
    let _frame = 0;

	const _update = (time) => {
		if( _mesh.material == null ) return;
        _frame ++;

        _mesh.material.uniforms.time.value = time * .001;
        _mesh.material.uniforms.frame.value = _frame;
        _mesh.material.uniforms.randomSeed.value = Math.random();
        _mesh.material.uniforms.fov.value = _camera.fov * Math.PI / 180;
       	_mesh.material.uniforms.raymarchMaximumDistance.value = _distance;
       	_mesh.material.uniforms.raymarchPrecision.value = _precision;
        _mesh.material.uniforms.camera.value = _camera.position;
        _mesh.material.uniforms.target.value = _target;

        _camera.lookAt( _target );
	}

	const _loadFragmentShader = async( fragmentUrl ) => {
		return promise = new Promise(resolve => {
			var req = new XMLHttpRequest();
	        req.open( "GET", fragmentUrl );
	        req.onload = function (e) {
	            return resolve( e.target.responseText );
	        };

	        req.send();
		})
    }

	const _setFragmentShader = (fragment) => {
		_startTime = Date.now();
        _mesh.material = new THREE.ShaderMaterial({

            uniforms :{
                resolution:{ type:"v2", value:new THREE.Vector2( this.width, this.height ) },
                time:{ type:"f", value:0 },
                frame:{ type:"i", value:0 },
                randomSeed:{ type:"f", value:Math.random() },
                fov:{ type:"f", value:45 },
                camera:{ type:"v3", value:_camera.position },
                target:{ type:"v3", value:_target },
                raymarchMaximumDistance:{ type:"f", value:_distance },
                raymarchPrecision:{ type:"f", value:_precision}

            },
            vertexShader : "void main() {gl_Position =  vec4( position, 1.0 );}",
            fragmentShader : fragment,
            transparent:true
        });
        
        _update();
	}

	return {
		
		onStart:async({distance, precision, fragmentUrl}) => {
			_distance = distance || 50;
			_precision = precision || 0.01;
			
			document.body.appendChild( _renderer.domElement );
			_width = window.innerWidth;
			_height =  window.innerHeight;
			_renderer.setSize( window.innerWidth, window.innerHeight );
        	_geom.setAttribute( 'position', new THREE.BufferAttribute( new Float32Array([   -1,-1,0, 1,-1,0, 1,1,0, -1, -1, 0, 1, 1, 0, -1, 1, 0]), 3 ) );
        	_scene.add( _mesh );
            //var helper = new THREE.PolarGridHelper(4, 8, 8, 48, new THREE.Color(0xaaaaaa), new THREE.Color(0xcccccc));
           // _scene.add(helper);
            
            _scene.background = new THREE.Color(0x000000);
        	const fragment = await _loadFragmentShader(fragmentUrl);

        	_setFragmentShader(fragment);

		},
	    setTexture:(name, url) => {
        	_mesh.material.uniforms[ name ] = {type:'t', value:null };
        	_tl.load( url, function(texture){

	            _mesh.material.uniforms[ name ].value = texture;
	            _mesh.material.needsUpdate = true;
	            
	            texture.needsUpdate = true;

	        });
	    },
	    setCubemap:( name, urls ) => {
        	_mesh.material.uniforms[ name ] = {type:'t', value:null };
        	_cl.load( urls, function(texture) {
            	_mesh.material.uniforms[ name ].value = texture;
            	_mesh.material.needsUpdate = true;
            
            	texture.needsUpdate = true;
        	});
    	},
    	setUniform:( name, type, value ) => {
    		if( _mesh.material == null ){
				throw new Error("raymarcher.getUniform: material not initialised, use setFragmentShader() first.");
				return null;
			}

        	_mesh.material.uniforms[ name ] = {type:type, value:value };
    	},
    	getUniform:( name ) =>{
			if( _mesh.material == null ){
				console.warn("raymarcher.getUniform: material not initialised, use setFragmentShader() first.");
				return null;
			}
        	return _mesh.material.uniforms[name];
    	},
    	setSize:( width, height ) => {

        	_width = width;
        	_height = height;

        	_renderer.setSize( width, height );

        	_camera.aspect = width / height;
        	_camera.updateProjectionMatrix();

	        if( _mesh.material != null ){
	            _mesh.material.uniforms.resolution.value.x = width;
	            _mesh.material.uniforms.resolution.value.y = height;
	        }
    	},
        addToScene:model => {
            _scene.add(model);
            model.visible=false;
        },
    	getCamera: _ => {
    		console.log(_camera);
    		return _camera;
    	},
    	getDomElement: _ => {
    		return _renderer.domElement;
    	},
    	getTarget:_ => {
    		return _target;
    	},
    	getMaterial:_ => {
    		return _mesh.material;
    	},
    	update:(time) => {
            _update(time);
            _renderer.render( _scene, _renderCamera );
	    }
	}
}