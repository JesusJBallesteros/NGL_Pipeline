function validationStr = fillValidators(propnames, props, namespacereg, className, inherited)
    validationStr = '';
    for i=1:length(propnames)
        nm = propnames{i};
        prop = props(nm);

        if (isa(prop, 'file.Attribute') || isa(prop, 'file.Dataset')) ...
                && prop.readonly && ~isempty(prop.value)
            % Need to add a validator for inherited and readonly properties. In 
            % the superclass these properties might not be read-only and due to
            % inheritance rules in MATLAB it is not possible to change property 
            % attributes of a property from public (in a superclass) to 
            % protected (in a subclass).
            if any(strcmp(nm, inherited))
                validationBody = fillReadOnlyValidator(nm, prop.value, className);
            else
                continue
            end
        elseif isa(prop, 'file.Link')
            validationBody = fillLinkValidation(nm, prop, namespacereg);
        else
            if startsWith(class(prop), 'file.')
                validationBody = fillUnitValidation(nm, prop, namespacereg);
            else % primitive type
                validationBody = fillDtypeValidation(nm, prop);
            end
        end

        headerStr = ['function val = validate_' nm '(obj, val)'];
        if isempty(validationBody)
            funcstionStr = [headerStr newline 'end'];
        else
            funcstionStr = strjoin({headerStr ...
                file.addSpaces(strtrim(validationBody), 4) 'end'}, newline);
        end
        validationStr = [validationStr newline funcstionStr];
    end
end

function unitValidationStr = fillUnitValidation(name, prop, namespaceReg)
    unitValidationStr = '';
    if ~isscalar(prop)
        constrained = cell(size(prop));
        for iProp = 1:length(prop)
            p = prop(iProp);
            assert(p.isConstrainedSet && ~isempty(p.type), 'Non-constrained anonymous set?');
            constrained{iProp} = ['''', namespaceReg.getFullClassName(p.type), ''''];
        end
        
        unitValidationStr = strjoin({...
            unitValidationStr, ...
            ['constrained = {', strjoin(constrained, ', '), '};'], ...
            ['types.util.checkSet(''', name, ''', struct(), constrained, val);'], ...
            }, newline);
    elseif isa(prop, 'file.Dataset')
        unitValidationStr = fillDatasetValidation(name, prop, namespaceReg);
    elseif isa(prop, 'file.Group')
        unitValidationStr = fillGroupValidation(name, prop, namespaceReg);
    elseif isa(prop, 'file.Attribute')
        unitValidationStr = strjoin({unitValidationStr...
            fillDtypeValidation(name, prop.dtype)...
            fillDimensionValidation(name, prop.shape)...
            }, newline);
    else % Link
        fullname = namespaceReg.getFullClassName(prop.type);
        unitValidationStr = fillDtypeValidation(name, fullname);
    end
end

function unitValidationStr = fillGroupValidation(name, prop, namespaceReg)
    if ~isempty(prop.type) && ~prop.isConstrainedSet
        fulltypename = namespaceReg.getFullClassName(prop.type);
        unitValidationStr = fillDtypeValidation(name, fulltypename);
        return;
    end

    namedprops = struct();
    constraints = {};
    if isempty(prop.type)
        %% process datasets
        % if type, check if constrained
        %   if constrained, add to constr
        %   otherwise, check type once
        % otherwise, check dtype
        for iDataset = 1:length(prop.datasets)
            dataset = p.datasets(iDataset);

            if isempty(dataset.type)
                namedprops.(dataset.name) = dataset.dtype;
            else
                type = namespaceReg.getFullClassName(dataset.type);
                if dataset.isConstrainedSet
                    constraints{end+1} = type;
                else
                    namedprops.(dataset.name) = type;
                end
            end
        end

        %% process groups
        % if type, check if constrained
        %   if constrained, add to constr
        %   otherwise, check type once
        % otherwise, error.  This shouldn't happen.
        for iSubGroup = 1:length(prop.subgroups)
            subGroup = prop.subgroups(iSubGroup);
            subGroupFullName = namespaceReg.getFullClassName(subGroup.type);
            assert(~isempty(subGroup.type), 'Weird case with two untyped groups');

            if isempty(subGroup.name)
                constraints{end+1} = subGroupFullName;
            else
                namedprops.(subGroup.name) = subGroupFullName;
            end
        end

        %% process attributes
        for iAttribute = 1:length(prop.attributes)
            Attribute = prop.attributes(iAttribute);
            namedprops.(Attribute.name) = Attribute.dtype;
        end

        %% process links
        for iLink = 1:length(prop.links)
            Link = prop.links(iLink);
            namespace = namespaceReg.getNamespace(Link.type);
            namedprops.(Link.name) = ['types.', namespace, '.', Link.type];
        end
    else
        constraints{end+1} = namespaceReg.getFullClassName(prop.type);
    end

    %% create unit validation string

    propnames = fieldnames(namedprops);
    unitValidationStr = 'namedprops = struct();';
    for i=1:length(propnames)
        nm = propnames{i};
        unitValidationStr = strjoin({unitValidationStr...
            ['namedprops.' nm ' = ''' namedprops.(nm) ''';']}, newline);
    end

    for iConstraint = 1:length(constraints)
        constraints{iConstraint} = ['''', constraints{iConstraint}, ''''];
    end

    unitValidationStr = strjoin({...
        unitValidationStr, ...
        ['constrained = {', strjoin(constraints, ','), '};'], ...
        ['types.util.checkSet(''', name, ''', namedprops, constrained, val);'], ...
        }, newline);
end

function unitValidationStr = fillDatasetValidation(name, prop, namespaceReg)
    unitValidationStr = '';
    if isempty(prop.type)
        unitValidationStr = strjoin({unitValidationStr...
            fillDtypeValidation(name, prop.dtype)...
            fillDimensionValidation(name, prop.shape)...
            }, newline);
    elseif prop.isConstrainedSet
        fullname = getFullClassName(namespaceReg, prop.type, name);
        if isempty(fullname)
            return
        end

        unitValidationStr = strjoin({unitValidationStr...
            ['constrained = { ''' fullname ''' };']...
            ['types.util.checkSet(''' name ''', struct(), constrained, val);']...
            }, newline);
    else
        fullname = getFullClassName(namespaceReg, prop.type, name);
        if isempty(fullname)
            return
        end
        unitValidationStr = [unitValidationStr newline fillDtypeValidation(name, fullname)];
    end
end

function validationStr = fillLinkValidation(name, prop, namespacereg)
    fullName = namespacereg.getFullClassName(prop.type);
    
    % Create a validation function body that 1) checks (validates) the
    % target if the input is a SoftLink type, otherwise 2) checks (validates) 
    % if the expected (target) type is provided. If the validation passes
    % and the value is not empty, it is wrapped in a SoftLink.

    validationStr = sprintf([ ...
        'if isa(val, ''types.untyped.SoftLink'')\n', ...
        '    if isprop(val, ''target'')\n', ...
        '        types.util.checkDtype(''%s'', ''%s'', val.target);\n', ...
        '    end\n', ...
        'else\n', ...
        '    %s\n', ...
        '    if ~isempty(val)\n', ...
        '        val = types.untyped.SoftLink(val);\n', ...
        '    end\n', ...
        'end' ...
        ], ...
        name, fullName, ...
        fillDtypeValidation(name, fullName) ...
    );
end

function fdvstr = fillDimensionValidation(name, shape)

    if iscell(shape)
        if ~isempty(shape) && iscell(shape{1})
            for i = 1:length(shape)
                for j = 1:length(shape{i})
                    shape{i}{j} = num2str(shape{i}{j});
                end
                shape{i} = ['[' strjoin(shape{i}, ',') ']'];
            end
            validShapeStr = ['{' strjoin(shape, ', ') '}'];
        else
            for i = 1:length(shape)
                shape{i} = num2str(shape{i});
            end
            validShapeStr = ['{[' strjoin(shape, ',') ']}'];
        end
    else
        validShapeStr = ['{[' num2str(shape) ']}'];
    end

    fdvstr = sprintf('types.util.validateShape(''%s'', %s, val)', name, validShapeStr);
end

%NOTE: can return empty strings
function fdvstr = fillDtypeValidation(name, type)
    if isstruct(type)
        fnames = fieldnames(type);
        fdvstr = strjoin({...
            'if isempty(val) || isa(val, ''types.untyped.DataStub'')'...
            '    return;'...
            'end'...
            'if ~istable(val) && ~isstruct(val) && ~isa(val, ''containers.Map'')'...
            ['    error(''NWB:Type:InvalidPropertyType'', ''Property `' name '` must be a table, struct, or containers.Map.'');']...
            'end'...
            'vprops = struct();'...
            }, newline);
        vprops = cell(length(fnames),1);
        for i=1:length(fnames)
            nm = fnames{i};
            if isa(type.(nm), 'containers.Map')
                %ref
                switch type.(nm)('reftype')
                    case 'region'
                        rt = 'RegionView';
                    case 'object'
                        rt = 'ObjectView';
                end
                typeval = ['types.untyped.' rt];
            else
                typeval = type.(nm);
            end
            vprops{i} = ['vprops.' nm ' = ''' typeval ''';'];
        end
        fdvstr = [fdvstr, newline, strjoin(vprops, newline), newline, ...
            'val = types.util.checkDtype(''' name ''', vprops, val);'];
    else
        fdvstr = '';
        if isa(type, 'containers.Map')
            %ref
            ref_t = type('reftype');
            switch ref_t
                case 'region'
                    rt = 'RegionView';
                case 'object'
                    rt = 'ObjectView';
            end
            ts = ['types.untyped.' rt];
            %there is no objective way to guarantee a reference refers to the
            %correct target type
            tt = type('target_type');
            fdvstr = ['% Reference to type `' tt '`' newline];
        elseif strcmp(type, 'any')
            fdvstr = '';
            return;
        else
            ts = strrep(type, '-', '_');
        end
        fdvstr = [fdvstr ...
            'val = types.util.checkDtype(''' name ''', ''' ts ''', val);'];
    end
end

function fdvstr = fillReadOnlyValidator(name, value, className)

    classNameSplit = strsplit(className, '.');
    shortName = classNameSplit{end};
    errorStr = sprintf( 'error(''NWB:Type:ReadOnlyProperty'', ''Unable to set the ''''%s'''' property of class ''''<a href="matlab:doc %s">%s</a>'''' because it is read-only.'')', name, className, shortName);  

    if ischar(value)
        condition = strjoin({ ...
            sprintf('if isequal(val, ''%s'')', value), ...
            sprintf('    val = ''%s'';', value ), ...
                    'else' }, newline);
    elseif isnumeric(value) || islogical(value)
        condition = strjoin({ ...
            sprintf('if isequal(val, %d)', value), ...
            sprintf('    val = %d;', value ), ...
                    'else' }, newline);
    else
        % Note: According to the documentation for Attribute specification keys
        % (https://schema-language.readthedocs.io/en/latest/description.html#sec-attributes-spec),
        % the above cases should be sufficient.
        error('NWB:ClassGenerator:ReadOnlyValidatorNotImplemented', ...
            'Read-only validator is not implemented for values of type "%s"', class(value))
    end
    
    fdvstr = strjoin({...
            condition, ...
            sprintf('    %s', errorStr), ...
            'end' }, newline );
end

function fullname = getFullClassName(namespaceReg, propType, name)
    fullname = '';
    try
        fullname = namespaceReg.getFullClassName(propType);
    catch ME
        if ~endsWith(ME.identifier, 'Namespace:NotFound')
            rethrow(ME);
        end

        warning('NWB:Fill:Validators:NamespaceNotFound',...
            ['Namespace could not be found for type `%s`.' ...
            '  Skipping Validation for property `%s`.'], propType, name);
    end
end
